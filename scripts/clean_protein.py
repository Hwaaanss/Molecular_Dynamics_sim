#!/usr/bin/env python3
# =============================================================================
# clean_protein.py --in IN.pdb --out OUT.pdb
#   PDBFixer 로 단백질 PDB 를 정제 + 결손 구조 복원.
#
#   결정구조에는 흔히 다음이 빠져 있어 pdb2gmx 가 실패하거나
#   비정상 결합을 만든다. 이를 pdb2gmx 이전에 보정한다:
#     - 결손된 무거운 원자(곁사슬)         → addMissingAtoms
#     - 결손된 잔기(내부 루프)             → addMissingResidues
#     - 비표준 잔기(MSE 등)               → replaceNonstandardResidues
#     - 물/공결정 리간드/이온(헤테로원자)   → removeHeterogens
#
#   수소는 추가하지 않는다 (pdb2gmx -ignh 가 포스필드 규칙대로 생성).
# =============================================================================
import argparse
import sys

from pdbfixer import PDBFixer
from openmm.app import PDBFile


def main():
    ap = argparse.ArgumentParser(description="PDBFixer 단백질 정제/결손복원")
    ap.add_argument("--in", dest="inp", required=True, help="입력 PDB")
    ap.add_argument("--out", required=True, help="출력(정제) PDB")
    ap.add_argument("--keep-water", action="store_true",
                    help="결정수 유지 (기본: 제거)")
    ap.add_argument("--fill-termini", action="store_true",
                    help="말단 결손 잔기까지 모델링 (기본: 내부 갭만 채움)")
    ap.add_argument("--no-minimize", action="store_true",
                    help="OpenMM 사전 최소화 생략")
    ap.add_argument("--minimize-iters", dest="minimize_iters", type=int, default=500,
                    help="OpenMM 최소화 반복 횟수 (기본 500)")
    args = ap.parse_args()

    print(f"[clean_protein] 입력: {args.inp}")
    fixer = PDBFixer(filename=args.inp)

    # 1) 물/리간드/이온 등 헤테로원자 제거 (기존 결합 리간드 제거 포함)
    fixer.removeHeterogens(keepWater=args.keep_water)

    # 2) 결손 잔기 탐색
    fixer.findMissingResidues()

    # 기본: 사슬 말단(꼬리) 결손은 모델링하지 않고 내부 갭만 채운다.
    #       (미해석 말단을 길게 모델링하면 부정확한 꼬리가 생기기 때문)
    if not args.fill_termini and fixer.missingResidues:
        chains = list(fixer.topology.chains())
        for key in list(fixer.missingResidues.keys()):
            chain_idx, res_pos = key
            n_res = len(list(chains[chain_idx].residues()))
            if res_pos == 0 or res_pos == n_res:   # 사슬 시작/끝
                del fixer.missingResidues[key]

    if fixer.missingResidues:
        print(f"[clean_protein] 내부 결손 잔기 복원: {fixer.missingResidues}")

    # 3) 비표준 잔기 치환 (예: MSE -> MET)
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        print(f"[clean_protein] 비표준 잔기 치환: {fixer.nonstandardResidues}")
    fixer.replaceNonstandardResidues()

    # 4) 결손 무거운 원자(곁사슬 등) 추가
    fixer.findMissingAtoms()
    if fixer.missingAtoms:
        print(f"[clean_protein] 결손 원자 보정: {len(fixer.missingAtoms)}개 잔기")
    fixer.addMissingAtoms()

    topology, positions = fixer.topology, fixer.positions

    # 5) OpenMM L-BFGS 사전 최소화 (선택)
    #    모델링된 루프/곁사슬 경계의 충돌은 GROMACS steepest descent 로 잘 안 풀려
    #    NVT 초반에 시스템이 터질 수 있다. OpenMM 최소화로 미리 완화한다.
    #    (실패하면 건너뛰고 미최소화 구조를 저장)
    if not args.no_minimize:
        try:
            topology, positions = openmm_minimize(fixer, args.minimize_iters)
            print(f"[clean_protein] OpenMM 사전 최소화 완료 ({args.minimize_iters} iters)")
        except Exception as e:  # noqa: BLE001
            print(f"[clean_protein] [warn] OpenMM 최소화 건너뜀: {e}")

    # 6) 저장 (수소 포함; pdb2gmx -ignh 가 재생성하므로 무방)
    with open(args.out, "w") as fh:
        PDBFile.writeFile(topology, positions, fh, keepIds=True)

    n_atom = sum(1 for _ in topology.atoms())
    n_res = sum(1 for _ in topology.residues())
    print(f"[clean_protein] 완료 → {args.out} ({n_res} residues, {n_atom} atoms)")


def openmm_minimize(fixer, iters):
    """PDBFixer 결과에 수소를 추가하고 진공에서 짧게 에너지 최소화.
       모델링된 영역의 원자 충돌을 완화한다. (topology, positions) 반환."""
    from openmm.app import ForceField, Modeller, NoCutoff, HBonds, Simulation
    from openmm import LangevinMiddleIntegrator, unit

    fixer.addMissingHydrogens(7.0)
    modeller = Modeller(fixer.topology, fixer.positions)
    ff = ForceField("amber14-all.xml")
    system = ff.createSystem(modeller.topology, nonbondedMethod=NoCutoff,
                             constraints=HBonds)
    integrator = LangevinMiddleIntegrator(300 * unit.kelvin,
                                          1.0 / unit.picosecond,
                                          0.001 * unit.picoseconds)
    sim = Simulation(modeller.topology, system, integrator)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=iters)
    state = sim.context.getState(getPositions=True)
    return modeller.topology, state.getPositions()


if __name__ == "__main__":
    main()
