# 체크포인트 / 재시작 가이드 (-cpi)

로컬 워크스테이션 특성상 정전·재부팅·강제종료로 시뮬레이션이 중단될 수 있다.
GROMACS 의 **체크포인트(.cpt)** 메커니즘으로 중단 지점부터 이어서 실행한다.

---

## 1. 체크포인트는 어떻게 생기나

`gmx mdrun` 은 기본적으로 **15분마다** `*.cpt` (예: `md.cpt`) 와 직전 백업
`*_prev.cpt` 를 자동 저장한다. 간격을 바꾸려면:

```bash
gmx mdrun -deffnm md -cpt 5 ...   # 5분마다 체크포인트
```

---

## 2. 중단된 MD 이어서 실행

핵심 옵션은 `-cpi <체크포인트>` 와 `-append`(기본 동작) 이다.

```bash
cd output/6VQN
gmx mdrun -deffnm md -cpi md.cpt \
  -nb gpu -pme gpu -bonded gpu -update gpu -ntmpi 1 -ntomp 8 -gpu_id 0
```

- `-cpi md.cpt` : 해당 체크포인트의 좌표·속도·스텝에서 재개
- `-append`     : 기존 `md.xtc`, `md.log`, `md.edr` 에 **이어쓰기** (기본값)
- 동일 `-deffnm` 을 그대로 사용해야 파일명이 일치한다

> 체크포인트 파일이 손상됐다면 직전 백업으로 재개:
> `gmx mdrun -deffnm md -cpi md_prev.cpt ...`

---

## 3. 이 저장소 스크립트의 자동 재시작

`scripts/03_run_md.sh` 의 `run_mdrun()` 래퍼가 다음을 자동 처리한다:

1. `<stage>.gro` 가 있으면 → **완료된 단계로 보고 건너뜀**
2. `<stage>.cpt` 가 있고 `.gro` 가 없으면 → **`-cpi` 로 자동 재시작**

따라서 중단 후에는 동일 명령을 다시 실행하기만 하면 된다:

```bash
./scripts/run_pipeline.sh 6VQN md      # 끊긴 지점부터 자동 재개
./scripts/run_pipeline.sh 6VQN         # 전체 재실행해도 완료분은 skip
```

---

## 4. 장시간 실행을 위한 백그라운드 구동

SSH 세션이 끊겨도 계속 돌도록 `nohup` 으로 실행하고 로그를 남긴다:

```bash
cd /home/dioxide421/Molecular_Dynamics_sim
nohup ./scripts/run_pipeline.sh 6VQN md > 6VQN_md.out 2>&1 &
echo $! > 6VQN_md.pid          # PID 기록 (필요시 kill)

# 진행 모니터링
tail -f 6VQN_md.out
tail -f output/6VQN/md.log
```

(가능하면 `tmux` / `screen` 세션 안에서 실행하는 것을 권장.)

---

## 5. 더 긴 시뮬레이션으로 연장 (.tpr 연장)

이미 끝난 production 을 더 늘리고 싶다면 `convert-tpr` 로 스텝을 연장한 뒤
`-cpi` 로 이어서 돌린다:

```bash
cd output/6VQN
gmx convert-tpr -s md.tpr -extend 2000 -o md_ext.tpr   # +2000 ps 연장
gmx mdrun -s md_ext.tpr -deffnm md -cpi md.cpt \
  -nb gpu -pme gpu -bonded gpu -update gpu -ntmpi 1 -ntomp 8 -gpu_id 0
```
