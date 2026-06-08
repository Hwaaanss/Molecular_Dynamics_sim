# Molecular_Dynamics_sim

하나의 약물(SDF)이 두 단백질 타겟(**6VQN**, **1VJY**)과 각각 결합한 복합체에 대해
GROMACS 분자동역학(MD) 시뮬레이션을 수행하는 자동화 파이프라인.
단백질 이름만 인자로 바꾸면 동일한 과정이 그대로 반복된다.

## 디렉터리 구조
```
Molecular_Dynamics_sim/
├── environment.yml              # conda 환경 'moledyn'
├── input/
│   ├── proteins/                # 6VQN.pdb, 1VJY.pdb
│   └── ligand/                  # drug.sdf
├── scripts/
│   ├── common.sh                # 공통 설정/헬퍼 (경로·포스필드·GPU)
│   ├── fetch_inputs.sh          # RCSB 에서 PDB 다운로드
│   ├── 00_prepare_ligand.py     # SDF → ACPYPE → 리간드 토폴로지
│   ├── 01_prepare_protein.sh    # PDB 정제 + pdb2gmx
│   ├── 02_build_system.sh       # 복합체 병합·박스·솔베이션·이온화
│   ├── 03_run_md.sh             # EM → NVT → NPT → Production (+체크포인트)
│   ├── 04_analyze.sh            # RMSD/RMSF → XVG
│   ├── 05_plot.py               # XVG → PNG
│   └── run_pipeline.sh          # 통합 드라이버 (단백질 이름 인자)
├── docs/
│   ├── GPU_SETUP.md             # A100 GPU 가속 빌드/실행
│   └── CHECKPOINT_RESTART.md    # 정전·중단 후 -cpi 재시작
└── output/<PROTEIN>/            # 단백질별 결과
```

## 사용법
### setting
```bash
conda create -n moledyn python=3.10 -y
tmux attach -t md
conda activate moledyn
```

위 최소 설정 대신, 필요한 도구를 한 번에 설치하려면 `environment.yml` 을 사용한다.
```bash
# 방법 A) environment.yml 로 한 번에 (권장)
conda env create -f environment.yml
conda activate moledyn

# 방법 B) 수동 설치
conda install -c conda-forge -c bioconda \
  rdkit openbabel ambertools acpype mdtraj gromacs matplotlib numpy pandas wget -y
```
> A100 GPU 가속을 최대로 쓰려면 GROMACS 를 CUDA 로 소스 빌드해야 한다.
> 절차와 환경변수는 [docs/GPU_SETUP.md](docs/GPU_SETUP.md) 참고.

### 1. 입력 파일 준비
```bash
# 단백질 PDB 자동 다운로드 (기본: 6VQN, 1VJY)
./scripts/fetch_inputs.sh
# 약물은 직접 배치
cp /path/to/your_drug.sdf input/ligand/drug.sdf
```

### 2. 전체 파이프라인 실행
단백질 이름만 바꿔 두 시스템을 각각 실행한다. (리간드 토폴로지는 첫 실행 시
1회 생성되어 두 시스템이 공유한다.)
```bash
./scripts/run_pipeline.sh 6VQN     # [6VQN + 약물] 전체 파이프라인
./scripts/run_pipeline.sh 1VJY     # [1VJY + 약물] 전체 파이프라인
```

### 3. 단계별 실행 (선택)
필요한 단계만 따로 실행할 수 있다.
```bash
./scripts/run_pipeline.sh 6VQN ligand    # 리간드 토폴로지만 (공유)
./scripts/run_pipeline.sh 6VQN protein   # 단백질 준비
./scripts/run_pipeline.sh 6VQN build     # 복합체·박스·솔베이션·이온
./scripts/run_pipeline.sh 6VQN md        # EM/NVT/NPT/Production
./scripts/run_pipeline.sh 6VQN analyze   # RMSD/RMSF + PNG
```
MD 단계는 4단계로 진행되며 각 단계의 `.mdp` 는 스크립트 내에서 자동 생성된다:
1. **Energy Minimization** (steep, emtol 1000)
2. **NVT 평형화** (100 ps, V-rescale 300 K, 위치제한)
3. **NPT 평형화** (100 ps, C-rescale 1 bar, 위치제한)
4. **Production MD** (테스트용 기본 1 ns; `PROD_NSTEPS` 로 조절)

### 4. 결과 확인
```
output/6VQN/analysis/
├── rmsd_backbone.png    # 단백질 backbone RMSD
├── rmsd_ligand.png      # 리간드 RMSD (단백질에 fit)
└── rmsf_backbone.png    # 잔기별 RMSF
```

### 5. 주요 환경변수 (override 가능)
| 변수 | 기본값 | 설명 |
|------|--------|------|
| `LIGAND_SDF` | `input/ligand/drug.sdf` | 약물 SDF 경로 |
| `FORCEFIELD` | `amber99sb-ildn` | pdb2gmx 포스필드 |
| `WATER_MODEL` | `tip3p` | 물 모델 |
| `LIG_ATOMTYPE` | `gaff2` | 리간드 atom type (gaff/gaff2) |
| `PROD_NSTEPS` | `500000` | production 스텝 수 (2 fs 기준 1 ns) |
| `GPU_ID` / `NT_OMP` / `NT_MPI` | `0` / `8` / `1` | GPU·스레드 자원 |

예) 2 ns production 을 GPU 1번에서:
```bash
PROD_NSTEPS=1000000 GPU_ID=1 ./scripts/run_pipeline.sh 6VQN md
```

### 6. 중단 후 재시작 / GPU 가속
- 정전·멈춤 후 동일 명령을 다시 실행하면 완료 단계는 건너뛰고 `.cpt`
  체크포인트에서 자동 재개된다. → [docs/CHECKPOINT_RESTART.md](docs/CHECKPOINT_RESTART.md)
- GPU 빌드/실행 환경변수 → [docs/GPU_SETUP.md](docs/GPU_SETUP.md)

### 7. 백그라운드 장시간 실행 예시
세션이 끊겨도 계속 돌도록 `tmux` 세션 안에서 실행하는 것을 권장한다.
```bash
tmux new -s md
nohup ./scripts/run_pipeline.sh 6VQN md > 6VQN_md.out 2>&1 &
tail -f output/6VQN/md.log     # 진행 모니터링 (Ctrl+b d 로 detach)
```
