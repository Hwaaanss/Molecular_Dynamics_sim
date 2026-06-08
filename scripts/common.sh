#!/usr/bin/env bash
# =============================================================================
# common.sh
#   Protein-Ligand MD 파이프라인 공통 설정 / 헬퍼 모음.
#   모든 스테이지 스크립트(01~05)와 run_pipeline.sh 에서 `source` 한다.
#
#   * 절대 단독 실행용이 아니다 (라이브러리).
#   * 경로/포스필드/GPU 자원은 환경변수로 override 가능하다.
# =============================================================================

# ----- bash 안전 옵션 --------------------------------------------------------
set -euo pipefail

# ----- 로깅 (이후 블록에서 사용하므로 최상단에 정의) -------------------------
log()  { printf '\033[1;34m[%s]\033[0m %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
warn() { printf '\033[1;33m[WARN]\033[0m %s\n' "$*" >&2; }
die()  { printf '\033[1;31m[ERROR]\033[0m %s\n' "$*" >&2; exit 1; }

# ----- 경로 ------------------------------------------------------------------
# 이 파일의 위치를 기준으로 프로젝트 루트를 계산한다.
COMMON_SH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${COMMON_SH_DIR}/.." && pwd)"

INPUT_DIR="${INPUT_DIR:-${ROOT_DIR}/input}"
OUTPUT_DIR="${OUTPUT_DIR:-${ROOT_DIR}/output}"
PROTEIN_INPUT_DIR="${PROTEIN_INPUT_DIR:-${INPUT_DIR}/proteins}"
LIGAND_INPUT_DIR="${LIGAND_INPUT_DIR:-${INPUT_DIR}/ligand}"

# 약물(SDF) 입력 파일. 두 단백질 시스템이 동일한 약물을 공유한다.
LIGAND_SDF="${LIGAND_SDF:-${LIGAND_INPUT_DIR}/drug.sdf}"
# ACPYPE 결과(공유)가 저장될 위치
LIGAND_WORK_DIR="${LIGAND_WORK_DIR:-${OUTPUT_DIR}/ligand}"
# 리간드 residue / moleculetype 이름
LIG_NAME="${LIG_NAME:-LIG}"

# ----- 시뮬레이션 파라미터 ---------------------------------------------------
FORCEFIELD="${FORCEFIELD:-amber99sb-ildn}"   # pdb2gmx -ff
WATER_MODEL="${WATER_MODEL:-tip3p}"          # pdb2gmx -water
BOX_TYPE="${BOX_TYPE:-dodecahedron}"         # editconf -bt
BOX_DIST="${BOX_DIST:-1.0}"                  # editconf -d (nm, 용질-박스 최소거리)
LIG_ATOMTYPE="${LIG_ATOMTYPE:-gaff2}"        # acpype -a (gaff / gaff2)
LIG_CHARGE_METHOD="${LIG_CHARGE_METHOD:-bcc}" # acpype -c (AM1-BCC)
LIG_PH="${LIG_PH:-7.4}"                       # openbabel 양성자화 pH
SALT_CONC="${SALT_CONC:-0.15}"                # genion -conc (mol/L)

# 평형화/Production 길이 (2 fs/step). 테스트용 기본값. 스모크테스트 시 줄여서 override.
NVT_NSTEPS="${NVT_NSTEPS:-50000}"     # 100 ps
NPT_NSTEPS="${NPT_NSTEPS:-50000}"     # 100 ps
PROD_NSTEPS="${PROD_NSTEPS:-500000}"  # 1 ns

# ----- GPU / 병렬 자원 (A100 80GB 기준 기본값) -------------------------------
GPU_ID="${GPU_ID:-0}"
NT_MPI="${NT_MPI:-1}"
NT_OMP="${NT_OMP:-8}"

# GROMACS 실행 바이너리.
#  CUDA 소스빌드(A100 가속)가 있으면 우선 사용하고, 없으면 PATH 의 gmx 사용.
#  (CUDA 빌드 절차는 docs/GPU_SETUP.md. 설치 위치는 GMX_CUDA_PREFIX 로 override 가능)
GMX_CUDA_PREFIX="${GMX_CUDA_PREFIX:-${HOME}/gromacs-gpu}"
if [[ -z "${GMX:-}" ]]; then
  if [[ -x "${GMX_CUDA_PREFIX}/bin/gmx" ]]; then
    GMX="${GMX_CUDA_PREFIX}/bin/gmx"
  else
    GMX="gmx"
  fi
fi

# ----- mdrun GPU 오프로딩 플래그 (백엔드 자동 감지) --------------------------
# 빌드에 따라 가능한 오프로딩이 다르다:
#   CUDA/SYCL : -nb -pme -bonded -update 모두 GPU 가능
#   OpenCL    : -nb, -pme 만 GPU (-bonded / -update 미지원)
#   (없음)    : CPU 전용
# `gmx --version` 의 "GPU support:" 항목으로 판별. GPU_BACKEND 로 override 가능.
GPU_BACKEND="${GPU_BACKEND:-$(${GMX} --version 2>/dev/null | awk -F: '/GPU support/{gsub(/[ \t]/,"",$2); print toupper($2)}')}"

case "${GPU_BACKEND}" in
  *CUDA*|*SYCL*)
    # 이 노드(드라이버 535 + 멀티 A100)에서 검증된 사실:
    #   -pme gpu / -update gpu → hang,  -bonded gpu → illegal memory access.
    # 안정적인 유일한 조합은 '비결합(nb)만 GPU' 다. 단, -nb gpu 만 주면 GROMACS
    # 가 PME/update 를 자동(auto)으로 GPU 에 올려 hang 되므로 명시적으로 cpu 지정.
    MDRUN_GPU="-nb gpu -pme cpu -bonded cpu -ntmpi ${NT_MPI} -ntomp ${NT_OMP} -gpu_id ${GPU_ID}"
    MDRUN_GPU_UPDATE="-nb gpu -pme cpu -bonded cpu -update cpu -ntmpi ${NT_MPI} -ntomp ${NT_OMP} -gpu_id ${GPU_ID}"
    ;;
  *OPENCL*)
    # OpenCL 빌드: bonded/update 미지원이고, NVIDIA GPU 는 GROMACS-OpenCL 과
    # 호환되지 않는 경우가 많다("incompatible devices"). 따라서 GPU 를 강제하지
    # 않고 자동선택(-nb auto)에 맡긴다 → 호환 GPU 있으면 사용, 없으면 CPU.
    # A100 가속을 제대로 쓰려면 CUDA 빌드 필요 (docs/GPU_SETUP.md).
    warn "GROMACS 가 OpenCL 빌드입니다. NVIDIA A100 가속은 CUDA 빌드가 필요합니다(docs/GPU_SETUP.md)."
    warn "  현재는 GPU 자동선택(미호환 시 CPU 폴백)으로 실행합니다."
    MDRUN_GPU="-ntmpi ${NT_MPI} -ntomp ${NT_OMP}"
    MDRUN_GPU_UPDATE="${MDRUN_GPU}"
    ;;
  *)
    # GPU 미지원 빌드 → CPU 실행
    MDRUN_GPU="-ntmpi ${NT_MPI} -ntomp ${NT_OMP}"
    MDRUN_GPU_UPDATE="${MDRUN_GPU}"
    ;;
esac

# 명령 존재 확인
require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "필수 명령을 찾을 수 없습니다: $1"
}

# 이미 산출물이 있으면 스테이지 skip 하도록 돕는 헬퍼
# usage: skip_if_done <output_file> <stage_name>  -> 있으면 return 0(skip)
skip_if_done() {
  local out="$1" name="$2"
  if [[ -f "$out" ]]; then
    log "[$name] 산출물 존재 → 건너뜀 ($out)"
    return 0
  fi
  return 1
}
