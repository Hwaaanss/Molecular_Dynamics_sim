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

# Production MD 길이: 테스트용 기본 1 ns (2 fs * 500000 = 1 ns)
PROD_NSTEPS="${PROD_NSTEPS:-500000}"

# ----- GPU / 병렬 자원 (A100 80GB 기준 기본값) -------------------------------
GPU_ID="${GPU_ID:-0}"
NT_MPI="${NT_MPI:-1}"
NT_OMP="${NT_OMP:-8}"

# mdrun GPU 오프로딩 플래그
#  - EM(steep)은 -update gpu 미지원 → MDRUN_GPU 사용
#  - NVT/NPT/Production 은 update 까지 GPU 로 → MDRUN_GPU_UPDATE 사용
MDRUN_GPU="-nb gpu -pme gpu -bonded gpu -ntmpi ${NT_MPI} -ntomp ${NT_OMP} -gpu_id ${GPU_ID}"
MDRUN_GPU_UPDATE="${MDRUN_GPU} -update gpu"

# GROMACS 실행 바이너리 (필요시 gmx_mpi 등으로 override)
GMX="${GMX:-gmx}"

# ----- 로깅 ------------------------------------------------------------------
log()  { printf '\033[1;34m[%s]\033[0m %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
warn() { printf '\033[1;33m[WARN]\033[0m %s\n' "$*" >&2; }
die()  { printf '\033[1;31m[ERROR]\033[0m %s\n' "$*" >&2; exit 1; }

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
