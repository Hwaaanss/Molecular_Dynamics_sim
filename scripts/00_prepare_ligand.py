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

    # 2) OpenBabel: SDF -> mol2 (수소 추가 + pH 양성자화)
    mol2 = os.path.join(outdir, f"{args.name}.mol2")
    run(["obabel", sdf, "-O", mol2, "-p", str(args.ph), "--partialcharge", "gasteiger"])

    # 3) ACPYPE 실행 (outdir 에서). 결과는 <name>.acpype/ 폴더로 생성됨.
    #    -i 입력, -b basename, -c charge model, -a atom type, -n net charge
    run([
        "acpype",
        "-i", mol2,
        "-b", args.name,
        "-c", args.charge_method,
        "-a", args.atomtype,
        "-n", str(charge),
        "-o", "gmx",
    ], cwd=outdir)

    acpype_dir = os.path.join(outdir, f"{args.name}.acpype")
    expected = [f"{args.name}_GMX.gro", f"{args.name}_GMX.itp", f"{args.name}_GMX.top"]
    missing = [f for f in expected if not os.path.isfile(os.path.join(acpype_dir, f))]
    if missing:
        sys.exit(f"[error] ACPYPE 결과 누락: {missing} (확인: {acpype_dir})")

    print(f"\n[done] 리간드 토폴로지 생성 완료 → {acpype_dir}")
    print("       다음 단계(02_build_system.sh)가 이 폴더를 사용합니다.")


if __name__ == "__main__":
    main()
