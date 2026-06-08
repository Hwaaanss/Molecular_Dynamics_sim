#!/usr/bin/env bash
# =============================================================================
# fetch_inputs.sh [PDB_ID ...]
#   RCSB PDB 에서 단백질 구조를 내려받아 input/proteins/ 에 저장한다.
#   인자가 없으면 기본 타겟(6VQN, 1VJY)을 받는다.
#   (약물 SDF 는 사용자가 input/ligand/drug.sdf 로 직접 준비해야 한다.)
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

require_cmd wget
mkdir -p "${PROTEIN_INPUT_DIR}"

PDB_IDS=("$@")
[[ ${#PDB_IDS[@]} -eq 0 ]] && PDB_IDS=(6VQN 1VJY)

for id in "${PDB_IDS[@]}"; do
  out="${PROTEIN_INPUT_DIR}/${id}.pdb"
  if [[ -f "${out}" ]]; then
    log "${id}.pdb 이미 존재 → skip"
    continue
  fi
  log "다운로드: ${id} → ${out}"
  wget -q "https://files.rcsb.org/download/${id}.pdb" -O "${out}" \
    || die "다운로드 실패: ${id}"
done

log "단백질 입력 준비 완료. 약물은 ${LIGAND_SDF} 로 직접 배치하세요."
