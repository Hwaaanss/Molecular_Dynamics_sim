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

log "[${PROTEIN}] 1) PDB 정제: 물(HOH)/헤테로원자(HETATM)/리간드 제거"
# 표준 단백질 원자(ATOM) 와 체인 종결(TER) 레코드만 보존.
#  → 결정수, 공결정 리간드, 이온 등이 모두 제거된다.
#  (참고: MSE 등 비표준 잔기는 HETATM 이라 함께 제거됨. 필요시 별도 처리.)
grep -E '^(ATOM|TER)' "${SRC_PDB}" > protein_clean.pdb
echo "END" >> protein_clean.pdb

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
