#!/usr/bin/env bash
# =============================================================================
# 01_prepare_protein.sh <PROTEIN_NAME>
#   단백질 PDB 정제 + pdb2gmx 토폴로지 생성.
#
#   입력 : input/proteins/<PROTEIN_NAME>.pdb
#   출력 : output/<PROTEIN_NAME>/
#            protein_clean.pdb   (물/헤테로원자 제거)
#            protein.gro         (pdb2gmx 좌표)
#            topol.top           (단백질 토폴로지)
#            posre.itp           (단백질 위치제한)
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

PROTEIN="${1:?사용법: $0 <PROTEIN_NAME> (예: 6VQN)}"
WORKDIR="${OUTPUT_DIR}/${PROTEIN}"
SRC_PDB="${PROTEIN_INPUT_DIR}/${PROTEIN}.pdb"

require_cmd "${GMX}"
[[ -f "${SRC_PDB}" ]] || die "단백질 PDB 가 없습니다: ${SRC_PDB}"
mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

# 멱등성: 이미 완료돼 있으면 재생성하지 않는다(재실행 시 topol/좌표 불일치 방지).
#  다시 만들려면 output/<PROTEIN> 를 지우고 실행.
if [[ -f topol.top && -f protein.gro ]]; then
  log "[${PROTEIN}] 단백질 토폴로지 이미 존재 → skip"
  exit 0
fi

log "[${PROTEIN}] 1) PDB 정제 + 결손 구조 복원"
# 결정구조는 곁사슬 원자/내부 루프가 빠진 경우가 많다. 이를 보정하지 않으면
# pdb2gmx 가 (a) 결손 원자 fatal error 를 내거나 (b) 갭을 가로질러 비정상
# 장거리 결합을 만든다. PDBFixer 로 사전 복원한다 (없으면 단순 grep 정제).
if python3 -c "import pdbfixer" >/dev/null 2>&1; then
  python3 "${SCRIPT_DIR}/clean_protein.py" --in "${SRC_PDB}" --out protein_clean.pdb
else
  warn "PDBFixer 미설치 → 단순 정제(결손 원자/잔기 복원 불가)."
  warn "  권장: conda install -n moledyn -c conda-forge pdbfixer openmm"
  grep -E '^(ATOM|TER)' "${SRC_PDB}" > protein_clean.pdb
  echo "END" >> protein_clean.pdb
fi

n_atom=$(grep -c '^ATOM' protein_clean.pdb || true)
[[ "${n_atom}" -gt 0 ]] || die "정제 후 ATOM 레코드가 0개입니다. PDB 형식을 확인하세요."
log "[${PROTEIN}] 정제 완료 (ATOM 레코드 ${n_atom}개)"

log "[${PROTEIN}] 2) pdb2gmx: 포스필드=${FORCEFIELD}, water=${WATER_MODEL}"
# -ignh : 입력 수소를 무시하고 포스필드 규칙대로 재생성 (PDB 수소 명명 불일치 회피)
"${GMX}" pdb2gmx \
  -f protein_clean.pdb \
  -o protein.gro \
  -p topol.top \
  -i posre.itp \
  -ff "${FORCEFIELD}" \
  -water "${WATER_MODEL}" \
  -ignh

log "[${PROTEIN}] 단백질 준비 완료 → ${WORKDIR}/{protein.gro, topol.top}"
