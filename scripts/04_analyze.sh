#!/usr/bin/env bash
# =============================================================================
# 04_analyze.sh <PROTEIN_NAME>
#   Production 궤적에서 RMSD / RMSF 계산 → XVG 추출 → PNG 자동저장.
#
#   선행조건: 03_run_md.sh 완료 (md.tpr, md.xtc)
#   산출물 (output/<PROTEIN_NAME>/analysis/):
#     md_center.xtc           PBC 보정 + centering 궤적
#     rmsd_backbone.xvg/.png  단백질 backbone RMSD
#     rmsd_ligand.xvg/.png    리간드 RMSD (단백질에 fit)
#     rmsf_backbone.xvg/.png  잔기별 RMSF
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

PROTEIN="${1:?사용법: $0 <PROTEIN_NAME>}"
WORKDIR="${OUTPUT_DIR}/${PROTEIN}"
ANADIR="${WORKDIR}/analysis"

require_cmd "${GMX}"
[[ -f "${WORKDIR}/md.tpr" && -f "${WORKDIR}/md.xtc" ]] || die "md.tpr/md.xtc 없음. 먼저 03_run_md.sh 실행."
mkdir -p "${ANADIR}"
cd "${WORKDIR}"

# ---------------------------------------------------------------------------
log "[${PROTEIN}] 1) trjconv: PBC 보정 + 단백질 centering"
# 1차: 분자 단위 PBC 처리 + 단백질 중심정렬 (입력그룹: Protein, 출력그룹: System)
printf "Protein\nSystem\n" | "${GMX}" trjconv \
  -s md.tpr -f md.xtc -o "${ANADIR}/md_center.xtc" \
  -center -pbc mol -ur compact

# ---------------------------------------------------------------------------
# 분석은 md.tpr 의 기본 그룹(Backbone, ${LIG_NAME})을 그대로 사용한다.
# (resname 을 LIG 로 통일했으므로 'LIG' 그룹이 기본으로 존재)
log "[${PROTEIN}] 2) RMSD (backbone, 단백질에 fit)"
printf "Backbone\nBackbone\n" | "${GMX}" rms \
  -s md.tpr -f "${ANADIR}/md_center.xtc" \
  -o "${ANADIR}/rmsd_backbone.xvg" -tu ns

log "[${PROTEIN}] 3) RMSD (리간드, 단백질 backbone 에 fit)"
printf "Backbone\n${LIG_NAME}\n" | "${GMX}" rms \
  -s md.tpr -f "${ANADIR}/md_center.xtc" \
  -o "${ANADIR}/rmsd_ligand.xvg" -tu ns

log "[${PROTEIN}] 4) RMSF (잔기별, backbone)"
printf "Backbone\n" | "${GMX}" rmsf \
  -s md.tpr -f "${ANADIR}/md_center.xtc" \
  -o "${ANADIR}/rmsf_backbone.xvg" -res

# ---------------------------------------------------------------------------
log "[${PROTEIN}] 5) XVG → PNG 변환"
python3 "${SCRIPT_DIR}/05_plot.py" \
  --title "${PROTEIN}" \
  --outdir "${ANADIR}" \
  "${ANADIR}/rmsd_backbone.xvg" \
  "${ANADIR}/rmsd_ligand.xvg" \
  "${ANADIR}/rmsf_backbone.xvg"

log "[${PROTEIN}] 분석 완료 → ${ANADIR}/*.png"
