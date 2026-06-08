#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh <PROTEIN_NAME> [STAGE]
#   하나의 단백질 타겟에 대한 전체 MD 파이프라인을 단일 명령으로 수행한다.
#   재사용성: 단백질 이름만 인자로 바꾸면 6VQN / 1VJY 둘 다 동일하게 처리.
#
#   STAGE (선택, 기본 all):
#     ligand   - 리간드 토폴로지 생성 (두 시스템 공유, 1회면 충분)
#     protein  - 01 단백질 준비
#     build    - 02 시스템 구축
#     md       - 03 MD 실행 (EM/NVT/NPT/Production)
#     analyze  - 04 분석(RMSD/RMSF) + 05 PNG
#     all      - ligand → protein → build → md → analyze 전부
#
#   예시:
#     ./scripts/run_pipeline.sh 6VQN
#     ./scripts/run_pipeline.sh 1VJY
#     ./scripts/run_pipeline.sh 6VQN md      # MD 단계만 (재시작 포함)
#
#   * 정전/멈춤 후 재실행해도 완료된 산출물은 건너뛰고 .cpt 로 이어서 진행한다.
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

PROTEIN="${1:?사용법: $0 <PROTEIN_NAME> [ligand|protein|build|md|analyze|all]}"
STAGE="${2:-all}"

# ---------------------------------------------------------------------------
stage_ligand() {
  if [[ -d "${LIGAND_WORK_DIR}/${LIG_NAME}.acpype" ]]; then
    log "리간드 토폴로지 이미 존재 → skip (${LIGAND_WORK_DIR}/${LIG_NAME}.acpype)"
    return 0
  fi
  log "=== 리간드 준비 (공유) ==="
  python3 "${SCRIPT_DIR}/00_prepare_ligand.py" \
    --sdf "${LIGAND_SDF}" \
    --outdir "${LIGAND_WORK_DIR}" \
    --name "${LIG_NAME}" \
    --atomtype "${LIG_ATOMTYPE}" \
    --charge-method "${LIG_CHARGE_METHOD}" \
    --ph "${LIG_PH}"
}

stage_protein() { log "=== [${PROTEIN}] 단백질 준비 ===";  bash "${SCRIPT_DIR}/01_prepare_protein.sh" "${PROTEIN}"; }
stage_build()   { log "=== [${PROTEIN}] 시스템 구축 ===";  bash "${SCRIPT_DIR}/02_build_system.sh"   "${PROTEIN}"; }
stage_md()      { log "=== [${PROTEIN}] MD 실행 ===";      bash "${SCRIPT_DIR}/03_run_md.sh"        "${PROTEIN}"; }
stage_analyze() { log "=== [${PROTEIN}] 분석/시각화 ===";  bash "${SCRIPT_DIR}/04_analyze.sh"       "${PROTEIN}"; }

# ---------------------------------------------------------------------------
case "${STAGE}" in
  ligand)  stage_ligand ;;
  protein) stage_protein ;;
  build)   stage_build ;;
  md)      stage_md ;;
  analyze) stage_analyze ;;
  all)
    stage_ligand
    stage_protein
    stage_build
    stage_md
    stage_analyze
    ;;
  *) die "알 수 없는 STAGE: ${STAGE} (ligand|protein|build|md|analyze|all)" ;;
esac

log "[${PROTEIN}] 파이프라인 단계 '${STAGE}' 완료."
