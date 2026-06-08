#!/usr/bin/env python3
# =============================================================================
# 00_prepare_ligand.py
#   SDF 약물 파일을 GROMACS 용 리간드 토폴로지로 변환한다.
#
#   파이프라인:
#     SDF --(RDKit)--> 형식전하(net charge) 자동 산출
#     SDF --(OpenBabel)--> 3D + pH 양성자화된 mol2
#     mol2 --(ACPYPE / AmberTools)--> GAFF/GAFF2 + AM1-BCC 토폴로지
#
#   산출물 (OUTDIR/<NAME>.acpype/):
#     <NAME>_GMX.gro , <NAME>_GMX.itp , <NAME>_GMX.top , posre_<NAME>.itp
#
#   두 단백질 시스템(6VQN, 1VJY)이 공유하므로 1회만 실행하면 된다.
# =============================================================================
import argparse
import os
import shutil
import subprocess
import sys


def run(cmd, **kw):
    """명령을 echo 후 실행. 실패 시 즉시 종료."""
    print(f"[run] {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True, **kw)


def force_mol2_resname(mol2_path, resname):
    """mol2 의 residue(subst) 이름을 resname 으로 통일.
       OpenBabel 은 기본 'UNL1' 을 쓰는데, 그러면 GROMACS 잔기명이 LIG 가
       아니라 UNL 이 되어 이후 resname 기반 그룹선택이 실패한다."""
    out, section = [], None
    with open(mol2_path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("@<TRIPOS>"):
                section = s
                out.append(line)
                continue
            f = s.split()
            if section == "@<TRIPOS>ATOM" and len(f) >= 9:
                f[7] = resname   # subst_name 열
                out.append("{:>7} {:<8}{:>10.4f}{:>10.4f}{:>10.4f} {:<6}{:>3} {:<8}{:>10.4f}\n".format(
                    int(f[0]), f[1], float(f[2]), float(f[3]), float(f[4]),
                    f[5], int(f[6]), f[7], float(f[8])))
            elif section == "@<TRIPOS>SUBSTRUCTURE" and len(f) >= 2:
                f[1] = resname
                out.append(" " + " ".join(f) + "\n")
            else:
                out.append(line)
    with open(mol2_path, "w") as fh:
        fh.writelines(out)
    print(f"[info] mol2 residue 이름을 '{resname}' 로 통일")


def detect_net_charge(sdf_path):
    """RDKit 로 SDF 의 형식전하 합을 계산. 실패하면 None 반환."""
    try:
        from rdkit import Chem
    except ImportError:
        print("[warn] RDKit 미설치 → net charge 자동산출 불가, --charge 사용 권장")
        return None

    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False, sanitize=True)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        print("[warn] RDKit 가 SDF 를 읽지 못함 → net charge 자동산출 실패")
        return None
    charge = Chem.GetFormalCharge(mol)
    print(f"[info] RDKit 형식전하(net charge) = {charge}")
    return charge


def detect_is_3d(sdf_path):
    """SDF 가 실제 3D 좌표인지 판정. 2D/평면이면 False.
       판정 불가(RDKit 미설치 등) 시 None 반환."""
    try:
        from rdkit import Chem
    except ImportError:
        return None
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False, sanitize=True)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None or mol.GetNumConformers() == 0:
        return None
    conf = mol.GetConformer()
    if not conf.Is3D():
        return False
    # z 좌표 분산이 사실상 0 이면 평면(2D)으로 간주
    zs = [conf.GetAtomPosition(i).z for i in range(mol.GetNumAtoms())]
    z_spread = max(zs) - min(zs)
    return z_spread > 0.1


def main():
    ap = argparse.ArgumentParser(description="SDF -> GROMACS 리간드 토폴로지 (ACPYPE)")
    ap.add_argument("--sdf", required=True, help="입력 약물 SDF 경로")
    ap.add_argument("--outdir", required=True, help="결과 저장 디렉터리")
    ap.add_argument("--name", default="LIG", help="리간드 residue/molecule 이름 (기본 LIG)")
    ap.add_argument("--charge", type=int, default=None,
                    help="net charge 강제 지정 (미지정 시 RDKit 자동산출)")
    ap.add_argument("--atomtype", default="gaff2", choices=["gaff", "gaff2"],
                    help="ACPYPE atom type (기본 gaff2)")
    ap.add_argument("--charge-method", default="bcc",
                    help="ACPYPE charge method (기본 bcc = AM1-BCC)")
    ap.add_argument("--ph", type=float, default=7.4,
                    help="OpenBabel 양성자화 pH (기본 7.4)")
    ap.add_argument("--gen3d", default="auto", choices=["auto", "always", "never"],
                    help="3D 좌표 재생성 여부. "
                         "auto=RDKit 로 2D 판정 시에만 재생성(기본), "
                         "always=항상 재생성(결합 포즈 손실), "
                         "never=재생성 안함(에너지 최소화로 충돌만 완화)")
    ap.add_argument("--mini-steps", type=int, default=2500,
                    help="OpenBabel 에너지 최소화 스텝 수 (원자 충돌 완화, 기본 2500)")
    args = ap.parse_args()

    sdf = os.path.abspath(args.sdf)
    outdir = os.path.abspath(args.outdir)
    if not os.path.isfile(sdf):
        sys.exit(f"[error] SDF 파일이 없습니다: {sdf}")
    os.makedirs(outdir, exist_ok=True)

    # 필수 외부 명령 확인
    for exe in ("obabel", "acpype"):
        if shutil.which(exe) is None:
            sys.exit(f"[error] 필수 명령을 찾을 수 없습니다: {exe} (conda env 'moledyn' 활성화?)")

    # 1) net charge 결정
    charge = args.charge
    if charge is None:
        charge = detect_net_charge(sdf)
    if charge is None:
        print("[warn] net charge 를 결정하지 못해 0 으로 가정합니다. (--charge 로 지정 가능)")
        charge = 0

    # 2) OpenBabel: SDF -> mol2 (수소 추가 + pH 양성자화 + 기하 완화)
    #    원자 충돌(특히 추가된 수소끼리 < 0.5 Å)이 있으면 ACPYPE 가 중단되므로,
    #    양성자화 후 반드시 기하 정리(3D 재생성 또는 에너지 최소화)를 수행한다.
    mol2 = os.path.join(outdir, f"{args.name}.mol2")

    do_gen3d = args.gen3d == "always"
    if args.gen3d == "auto":
        is3d = detect_is_3d(sdf)
        if is3d is False:
            print("[info] 입력 SDF 가 2D/평면으로 판정 → 3D 좌표 재생성(--gen3d)")
            do_gen3d = True
        elif is3d is None:
            print("[warn] 3D 여부 판정 불가(RDKit) → 에너지 최소화만 수행")

    cmd = ["obabel", sdf, "-O", mol2, "-p", str(args.ph)]
    if do_gen3d:
        # 2D 또는 강제 재생성: 전체 3D 좌표 생성 (결합 포즈는 보존되지 않음)
        cmd += ["--gen3d"]
    else:
        # 3D 포즈 유지: 무거운 원자 위치를 출발점으로 짧은 에너지 최소화 → 충돌 해소
        cmd += ["--minimize", "--sd", "--steps", str(args.mini_steps), "--ff", "MMFF94"]
    run(cmd)

    # 잔기 이름을 LIG 로 통일 (이후 resname 기반 인덱스/분석이 일관되게 동작)
    force_mol2_resname(mol2, args.name)

    # 3) ACPYPE 실행 (outdir 에서). 결과는 <name>.acpype/ 폴더로 생성됨.
    #    -i 입력, -b basename, -c charge model, -a atom type, -n net charge
    try:
        run([
            "acpype",
            "-i", mol2,
            "-b", args.name,
            "-c", args.charge_method,
            "-a", args.atomtype,
            "-n", str(charge),
            "-o", "gmx",
        ], cwd=outdir)
    except subprocess.CalledProcessError:
        sys.exit(
            "[error] ACPYPE 실패.\n"
            "  'Atoms TOO close' 라면 리간드 기하에 원자 충돌이 남아 있는 것입니다.\n"
            "  해결: --gen3d always 로 3D 좌표를 새로 생성하거나(포즈 손실),\n"
            "        --mini-steps 값을 늘려 최소화를 강화하세요.\n"
            "  예) python3 scripts/00_prepare_ligand.py --sdf <SDF> --outdir output/ligand --gen3d always"
        )

    acpype_dir = os.path.join(outdir, f"{args.name}.acpype")
    expected = [f"{args.name}_GMX.gro", f"{args.name}_GMX.itp", f"{args.name}_GMX.top"]
    missing = [f for f in expected if not os.path.isfile(os.path.join(acpype_dir, f))]
    if missing:
        sys.exit(f"[error] ACPYPE 결과 누락: {missing} (확인: {acpype_dir})")

    print(f"\n[done] 리간드 토폴로지 생성 완료 → {acpype_dir}")
    print("       다음 단계(02_build_system.sh)가 이 폴더를 사용합니다.")


if __name__ == "__main__":
    main()
