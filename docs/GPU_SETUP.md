# GROMACS GPU 가속 설정 (NVIDIA A100 80GB)

이 문서는 GROMACS 가 A100 GPU 를 최대한 활용하도록 **빌드** 하고 **실행** 하는
방법을 정리한다.

---

## 1. 왜 소스 빌드인가

`conda` 의 `gromacs` 패키지는 대체로 CPU/일반 빌드라 GPU 오프로딩이 제한적이다.
전처리(`pdb2gmx`, `editconf`, `solvate`, `grompp`) 와 문법검증에는 충분하지만,
실제 production MD 를 A100 으로 가속하려면 **CUDA 지원 소스 빌드**를 권장한다.

---

## 2. 소스 빌드 (CUDA)

사전 요구: CUDA Toolkit, `cmake (>=3.18)`, `gcc`, MPI(선택).

```bash
# 0) CUDA 설치 확인
nvidia-smi
nvcc --version

# 1) 소스 내려받기 (버전은 환경에 맞게)
VER=2024.3
wget https://ftp.gromacs.org/gromacs/gromacs-${VER}.tar.gz
tar xfz gromacs-${VER}.tar.gz
cd gromacs-${VER}
mkdir build && cd build

# 2) CMake 구성: CUDA 가속 활성화
cmake .. \
  -DGMX_GPU=CUDA \
  -DGMX_CUDA_TARGET_SM=80 \          # A100 = compute capability 8.0
  -DCMAKE_INSTALL_PREFIX=$HOME/gromacs-gpu \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DGMX_SIMD=AVX2_256 \
  -DREGRESSIONTEST_DOWNLOAD=ON

# 3) 빌드 / (선택)검증 / 설치
make -j$(nproc)
make check          # 선택: 회귀테스트
make install

# 4) 환경 로드 (셸 세션마다 또는 ~/.bashrc 에 추가)
source $HOME/gromacs-gpu/bin/GMXRC
```

빌드 후 `gmx --version` 출력에서 `GPU support: CUDA` 를 확인한다.

---

## 3. 실행 시 환경변수

A100 단일 GPU + 다중 OpenMP 스레드 기준 권장 설정.

```bash
# 사용할 GPU 선택 (0번 GPU)
export CUDA_VISIBLE_DEVICES=0

# OpenMP 스레드 수 (mdrun -ntomp 와 일치시킨다)
export OMP_NUM_THREADS=8

# (선택) GPU 직접통신 최적화 — 단일 GPU 에서도 update/buffer-ops 가속
export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=true
```

이 저장소의 스크립트는 위 변수들을 override 할 수 있도록 작성돼 있다
(`scripts/common.sh` 의 `GPU_ID`, `NT_MPI`, `NT_OMP`).

---

## 4. mdrun GPU 오프로딩 플래그

`scripts/common.sh` 에서 자동 구성되는 플래그:

```bash
# NVT / NPT / Production (모든 항을 GPU 로)
gmx mdrun -deffnm md \
  -nb gpu -pme gpu -bonded gpu -update gpu \
  -ntmpi 1 -ntomp 8 -gpu_id 0

# Energy Minimization (steep 는 -update gpu 미지원 → 비결합력만)
gmx mdrun -deffnm em \
  -nb gpu -pme gpu -bonded gpu \
  -ntmpi 1 -ntomp 8 -gpu_id 0
```

각 항목 의미:
- `-nb gpu`     : 비결합 상호작용(non-bonded) GPU
- `-pme gpu`    : PME 정전기 GPU (단일 GPU 에 적합)
- `-bonded gpu` : 결합 상호작용 GPU
- `-update gpu` : 적분/제약(update+constraints) GPU → CPU-GPU 전송 최소화

> 환경변수 `GPU_ID`, `NT_OMP`, `NT_MPI` 로 자원을 조정할 수 있다.
> 예) `NT_OMP=16 GPU_ID=1 ./scripts/run_pipeline.sh 6VQN md`

---

## 5. 성능 점검

실행 후 `md.log` 끝부분의 **Performance (ns/day)** 와 GPU 활용도를 확인:

```bash
grep -A2 "Performance" output/6VQN/md.log
nvidia-smi dmon -s u    # 실시간 GPU 사용률
```

PME 가 병목이면 `-pme gpu` 유지한 채 `-ntomp` 를 늘려보고, GPU 메모리 부족 시
(드물지만 큰 박스에서) 박스 크기(`BOX_DIST`)를 점검한다.
