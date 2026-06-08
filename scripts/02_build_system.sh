#!/usr/bin/env bash
# =============================================================================
# 02_build_system.sh <PROTEIN_NAME>
#   단백질 + 리간드 복합체 구성 → 박스 → 솔베이션 → 이온화.
#
#   선행조건:
#     - 01_prepare_protein.sh 완료 (protein.gro, topol.top, posre.itp)
#     - 00_prepare_ligand.py 완료 (output/ligand/<LIG>.acpype/)
#
#   주요 산출물 (output/<PROTEIN_NAME>/):
#     complex.gro          단백질+리간드 병합 좌표
#     topol.top            리간드 #include + 분자수 반영
#     LIG.itp / atomtypes_LIG.itp / posre_LIG.itp
#     solv_ions.gro        솔베이션+이온화 완료 (EM 입력)
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

PROTEIN="${1:?사용법: $0 <PROTEIN_NAME>}"
WORKDIR="${OUTPUT_DIR}/${PROTEIN}"
ACPYPE_DIR="${LIGAND_WORK_DIR}/${LIG_NAME}.acpype"

require_cmd "${GMX}"
[[ -d "${WORKDIR}" ]]    || die "단백질 작업폴더가 없습니다. 먼저 01_prepare_protein.sh 실행: ${WORKDIR}"
[[ -f "${WORKDIR}/protein.gro" ]] || die "protein.gro 없음. 01_prepare_protein.sh 먼저 실행."
[[ -d "${ACPYPE_DIR}" ]] || die "리간드 ACPYPE 폴더 없음. 먼저 00_prepare_ligand.py 실행: ${ACPYPE_DIR}"
cd "${WORKDIR}"

LIG_GRO="${ACPYPE_DIR}/${LIG_NAME}_GMX.gro"
LIG_ITP_SRC="${ACPYPE_DIR}/${LIG_NAME}_GMX.itp"
LIG_TOP_SRC="${ACPYPE_DIR}/${LIG_NAME}_GMX.top"
POSRE_LIG_SRC="${ACPYPE_DIR}/posre_${LIG_NAME}.itp"

# ---------------------------------------------------------------------------
log "[${PROTEIN}] 1) 리간드 토폴로지 가공 (atomtypes 분리)"
# ACPYPE 버전에 따라 [ atomtypes ] 가 _GMX.itp 또는 _GMX.top 에 들어있다.
# GROMACS topol.top 에서는 atomtypes 가 어떤 moleculetype 보다 먼저,
# forcefield.itp include 직후에 와야 한다. 따라서:
#   (a) atomtypes 블록을 atomtypes_LIG.itp 로 추출
#   (b) LIG.itp 에서는 atomtypes 블록을 제거 (중복정의 방지)
if grep -q '^\[ *atomtypes *\]' "${LIG_ITP_SRC}"; then
  ATOMTYPE_SRC="${LIG_ITP_SRC}"
else
  ATOMTYPE_SRC="${LIG_TOP_SRC}"
fi

# (a) atomtypes 블록만 추출
awk '
  /^\[ *atomtypes *\]/        { grab=1 }
  grab && /^\[/ && !/atomtypes/ { grab=0 }
  grab                         { print }
' "${ATOMTYPE_SRC}" > atomtypes_LIG.itp
[[ -s atomtypes_LIG.itp ]] || die "리간드 [ atomtypes ] 추출 실패 (${ATOMTYPE_SRC})"

# (b) moleculetype itp 에서 atomtypes 제거
awk '
  /^\[ *atomtypes *\]/ { skip=1; next }
  /^\[/                { skip=0 }
  !skip                { print }
' "${LIG_ITP_SRC}" > LIG.itp

cp -f "${LIG_GRO}" LIG.gro
[[ -f "${POSRE_LIG_SRC}" ]] && cp -f "${POSRE_LIG_SRC}" "posre_${LIG_NAME}.itp"

# ACPYPE 는 posre include 를 _GMX.top 에만 넣고 .itp 에는 넣지 않는다.
# → LIG.itp 에 직접 추가해야 -DPOSRES_LIG 가 실제로 리간드를 제한한다.
if [[ -f "posre_${LIG_NAME}.itp" ]] && ! grep -q "posre_${LIG_NAME}.itp" LIG.itp; then
  cat >> LIG.itp << EOF

; Ligand position restraints
#ifdef POSRES_${LIG_NAME}
#include "posre_${LIG_NAME}.itp"
#endif
EOF
fi

# ---------------------------------------------------------------------------
log "[${PROTEIN}] 2) 단백질 + 리간드 좌표 병합 → complex.gro"
# .gro 구조: 1=title, 2=원자수, 3..N+2=원자, 마지막=박스벡터
np=$(sed -n '2p' protein.gro | tr -d '[:space:]')
nl=$(sed -n '2p' LIG.gro     | tr -d '[:space:]')
total=$(( np + nl ))
{
  echo "${PROTEIN} + ${LIG_NAME} complex"
  printf '%5d\n' "${total}"
  tail -n +3 protein.gro | head -n -1   # 단백질 원자 (박스줄 제외)
  tail -n +3 LIG.gro     | head -n -1   # 리간드 원자 (박스줄 제외)
  tail -n 1 protein.gro                 # 박스 벡터 (단백질 것 사용)
} > complex.gro
log "[${PROTEIN}] 병합 완료 (단백질 ${np} + 리간드 ${nl} = ${total} atoms)"

# ---------------------------------------------------------------------------
log "[${PROTEIN}] 3) topol.top 갱신 (#include + 분자수)"
# 멱등성: 이미 추가돼 있으면 다시 넣지 않는다.
if ! grep -q 'atomtypes_LIG.itp' topol.top; then
  # forcefield.itp include 직후에 리간드 atomtypes + itp include 삽입
  awk '
    { print }
    /#include.*forcefield\.itp/ && !done {
      print ""
      print "; Include ligand atom types"
      print "#include \"atomtypes_LIG.itp\""
      print "; Include ligand topology"
      print "#include \"LIG.itp\""
      done=1
    }
  ' topol.top > topol.top.tmp && mv topol.top.tmp topol.top
fi

# [ molecules ] 섹션에 리간드 1개 추가 (중복 방지)
if ! grep -qE "^\s*${LIG_NAME}\s+[0-9]+" topol.top; then
  printf '%-15s 1\n' "${LIG_NAME}" >> topol.top
fi

# ---------------------------------------------------------------------------
log "[${PROTEIN}] 4) editconf: ${BOX_TYPE} 박스 (-d ${BOX_DIST} nm)"
"${GMX}" editconf -f complex.gro -o box.gro -c -d "${BOX_DIST}" -bt "${BOX_TYPE}"

log "[${PROTEIN}] 5) solvate: 물 채우기 (${WATER_MODEL})"
"${GMX}" solvate -cp box.gro -cs spc216.gro -p topol.top -o solv.gro

# ---------------------------------------------------------------------------
log "[${PROTEIN}] 6) genion: 중성화 + ${SALT_CONC} M 염 추가"
# ions.mdp 즉석 생성 (genion 용 더미 tpr 생성에만 사용)
cat > ions.mdp << 'EOF'
; ions.mdp - genion 용 더미 파라미터
integrator  = steep
emtol       = 1000.0
nsteps      = 5000
nstlist     = 10
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
EOF

# maxwarn 2: 이온 추가 전이라 (1) net charge + Ewald 경고는 필연적이고,
#            (2) EM 전이라 '제외원자 거리>컷오프' 경고가 날 수 있다(EM 이 해소).
"${GMX}" grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2

# SOL 그룹을 이온으로 치환. 입력 그룹 선택은 stdin 으로 "SOL".
echo "SOL" | "${GMX}" genion \
  -s ions.tpr -o solv_ions.gro -p topol.top \
  -pname NA -nname CL -neutral -conc "${SALT_CONC}"

log "[${PROTEIN}] 시스템 구축 완료 → ${WORKDIR}/solv_ions.gro"
