#!/usr/bin/env bash
# =============================================================================
# 03_run_md.sh <PROTEIN_NAME>
#   에너지최소화(EM) → NVT → NPT → Production MD.
#   모든 .mdp 파일을 스크립트 내에서 즉석 생성(cat << EOF)한다.
#   GPU 오프로딩 및 체크포인트(-cpi) 재시작을 지원한다.
#
#   선행조건: 02_build_system.sh 완료 (solv_ions.gro, topol.top)
#   산출물 (output/<PROTEIN_NAME>/):
#     em.gro, nvt.gro, npt.gro, md.gro, md.xtc, md.tpr, *.cpt
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

PROTEIN="${1:?사용법: $0 <PROTEIN_NAME>}"
WORKDIR="${OUTPUT_DIR}/${PROTEIN}"

require_cmd "${GMX}"
[[ -f "${WORKDIR}/solv_ions.gro" ]] || die "solv_ions.gro 없음. 먼저 02_build_system.sh 실행."
cd "${WORKDIR}"

# ---------------------------------------------------------------------------
# 체크포인트 인지형 mdrun 래퍼.
#   - <deffnm>.gro 가 이미 있으면 완료된 것으로 보고 skip
#   - <deffnm>.cpt 가 있으면 -cpi 로 이어서 실행 (정전/멈춤 복구)
# 사용: run_mdrun <deffnm> <gpu_flags...>
# ---------------------------------------------------------------------------
run_mdrun() {
  local deffnm="$1"; shift
  local gpu_flags=("$@")
  if [[ -f "${deffnm}.gro" ]]; then
    log "[${PROTEIN}] ${deffnm} 이미 완료 → skip"
    return 0
  fi
  local cpi=()
  if [[ -f "${deffnm}.cpt" ]]; then
    log "[${PROTEIN}] ${deffnm}.cpt 발견 → 체크포인트에서 재시작"
    cpi=(-cpi "${deffnm}.cpt")
  fi
  "${GMX}" mdrun -deffnm "${deffnm}" "${cpi[@]}" "${gpu_flags[@]}"
}

# =========================== 1) 에너지 최소화 ===============================
log "[${PROTEIN}] === Stage 1: Energy Minimization ==="
cat > em.mdp << 'EOF'
; em.mdp - steepest descent 에너지 최소화
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000
nstlist         = 10
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
EOF

"${GMX}" grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1
# EM(steep)은 비동역학 적분기 → PME/update GPU 미지원, 비결합력만 GPU 오프로딩
run_mdrun em ${MDRUN_GPU}

# ---------------------------------------------------------------------------
# 온도커플링용 인덱스 생성: Protein_LIG / Water_and_ions
#   gmx select 로 그룹번호에 의존하지 않고 견고하게 생성.
# ---------------------------------------------------------------------------
if [[ ! -f index.ndx ]]; then
  log "[${PROTEIN}] index.ndx 생성 (Protein_LIG / Water_and_ions)"
  # 커스텀 인덱스를 주면 기본 그룹이 자동 포함되지 않으므로 System 도 명시한다
  # (production md.mdp 의 compressed-x-grps = System 에서 필요).
  "${GMX}" select -s em.tpr -on index.ndx \
    -select "\"System\" all" \
            "\"Protein_LIG\" group \"Protein\" or resname ${LIG_NAME}" \
            "\"Water_and_ions\" not (group \"Protein\" or resname ${LIG_NAME})"
fi

# =========================== 2) NVT 평형화 =================================
log "[${PROTEIN}] === Stage 2: NVT equilibration (${NVT_NSTEPS} steps) ==="
cat > nvt.mdp << EOF
; nvt.mdp - 위치제한 하 온도 평형화
define          = -DPOSRES -DPOSRES_${LIG_NAME}
integrator      = md
nsteps          = ${NVT_NSTEPS}   ; 2 fs/step
dt              = 0.002
; 출력
nstxout         = 0
nstvout         = 0
nstenergy       = 500
nstlog          = 500
; 결합 제약
continuation    = no
constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4
; 이웃탐색
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 20
rcoulomb        = 1.0
rvdw            = 1.0
; 정전기
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16
; 온도커플링
tcoupl          = V-rescale
tc-grps         = Protein_LIG Water_and_ions
tau_t           = 0.1 0.1
ref_t           = 300 300
; 압력커플링 없음
pcoupl          = no
; PBC / 분산보정
pbc             = xyz
DispCorr        = EnerPres
; 초기속도 생성
gen_vel         = yes
gen_temp        = 300
gen_seed        = -1
EOF

"${GMX}" grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 1
run_mdrun nvt ${MDRUN_GPU_UPDATE}

# =========================== 3) NPT 평형화 =================================
log "[${PROTEIN}] === Stage 3: NPT equilibration (${NPT_NSTEPS} steps) ==="
cat > npt.mdp << EOF
; npt.mdp - 위치제한 하 압력 평형화
define          = -DPOSRES -DPOSRES_${LIG_NAME}
integrator      = md
nsteps          = ${NPT_NSTEPS}   ; 2 fs/step
dt              = 0.002
; 출력
nstxout         = 0
nstvout         = 0
nstenergy       = 500
nstlog          = 500
; 결합 제약 (NVT 속도 이어받음)
continuation    = yes
constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4
; 이웃탐색
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 20
rcoulomb        = 1.0
rvdw            = 1.0
; 정전기
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16
; 온도커플링
tcoupl          = V-rescale
tc-grps         = Protein_LIG Water_and_ions
tau_t           = 0.1 0.1
ref_t           = 300 300
; 압력커플링 (C-rescale: 평형화에 적합)
pcoupl          = C-rescale
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = 1.0
compressibility = 4.5e-5
refcoord_scaling = com
; PBC / 분산보정
pbc             = xyz
DispCorr        = EnerPres
; 속도 생성 안함
gen_vel         = no
EOF

# NPT 는 NVT 의 체크포인트(속도)를 이어받는다 (-t nvt.cpt)
"${GMX}" grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn 1
run_mdrun npt ${MDRUN_GPU_UPDATE}

# =========================== 4) Production MD ==============================
prod_ns=$(awk "BEGIN{printf \"%.2f\", ${PROD_NSTEPS}*0.002/1000}")
log "[${PROTEIN}] === Stage 4: Production MD (${PROD_NSTEPS} steps ≈ ${prod_ns} ns) ==="
cat > md.mdp << EOF
; md.mdp - 위치제한 없는 production MD
integrator      = md
nsteps          = ${PROD_NSTEPS}
dt              = 0.002
; 출력 (압축 좌표 10 ps 간격)
nstxout         = 0
nstvout         = 0
nstenergy       = 5000
nstlog          = 5000
nstxout-compressed = 5000
compressed-x-grps  = System
; 결합 제약
continuation    = yes
constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4
; 이웃탐색
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 20
rcoulomb        = 1.0
rvdw            = 1.0
; 정전기
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16
; 온도커플링
tcoupl          = V-rescale
tc-grps         = Protein_LIG Water_and_ions
tau_t           = 0.1 0.1
ref_t           = 300 300
; 압력커플링 (Parrinello-Rahman: production 에 적합)
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = 1.0
compressibility = 4.5e-5
; PBC / 분산보정
pbc             = xyz
DispCorr        = EnerPres
; 속도 생성 안함
gen_vel         = no
EOF

"${GMX}" grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr -maxwarn 1
run_mdrun md ${MDRUN_GPU_UPDATE}

log "[${PROTEIN}] MD 완료 → ${WORKDIR}/md.xtc , md.gro , md.tpr"
