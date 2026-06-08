#!/usr/bin/env python3
# =============================================================================
# 05_plot.py --outdir DIR [--title T] file1.xvg [file2.xvg ...]
#   GROMACS XVG 파일을 파싱하여 PNG 그래프로 저장한다.
#   '@'(grace 명령) 와 '#'(주석) 라인을 건너뛰고, 축/제목 메타데이터를 활용한다.
# =============================================================================
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")  # headless 환경 (GUI 없음)
import matplotlib.pyplot as plt


def parse_xvg(path):
    """XVG -> (data[list of rows], meta dict). 숫자행만 data 에 수집."""
    meta = {"title": "", "xaxis": "x", "yaxis": "y", "legends": []}
    rows = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                continue
            if s.startswith("@"):
                # grace 메타데이터 파싱
                if "title" in s and "subtitle" not in s:
                    meta["title"] = _between_quotes(s) or meta["title"]
                elif "xaxis" in s and "label" in s:
                    meta["xaxis"] = _between_quotes(s) or meta["xaxis"]
                elif "yaxis" in s and "label" in s:
                    meta["yaxis"] = _between_quotes(s) or meta["yaxis"]
                elif s.startswith("@ s") and "legend" in s:
                    leg = _between_quotes(s)
                    if leg:
                        meta["legends"].append(leg)
                continue
            # 데이터 라인
            parts = s.split()
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue
    return rows, meta


def _between_quotes(s):
    a = s.find('"')
    b = s.rfind('"')
    return s[a + 1:b] if (a != -1 and b > a) else ""


def plot_one(path, outdir, title_prefix):
    rows, meta = parse_xvg(path)
    if not rows:
        print(f"[warn] 데이터 없음, 건너뜀: {path}")
        return None

    ncol = len(rows[0])
    xs = [r[0] for r in rows]

    fig, ax = plt.subplots(figsize=(8, 5))
    for c in range(1, ncol):
        ys = [r[c] for r in rows]
        label = meta["legends"][c - 1] if c - 1 < len(meta["legends"]) else f"col{c}"
        ax.plot(xs, ys, lw=1.2, label=label)

    ax.set_xlabel(meta["xaxis"])
    ax.set_ylabel(meta["yaxis"])
    base = os.path.splitext(os.path.basename(path))[0]
    title = f"{title_prefix} - {meta['title'] or base}" if title_prefix else (meta["title"] or base)
    ax.set_title(title)
    if ncol > 2 or meta["legends"]:
        ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    png = os.path.join(outdir, base + ".png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    print(f"[saved] {png}")
    return png


def main():
    ap = argparse.ArgumentParser(description="GROMACS XVG -> PNG 그래프")
    ap.add_argument("xvg", nargs="+", help="입력 XVG 파일들")
    ap.add_argument("--outdir", required=True, help="PNG 출력 디렉터리")
    ap.add_argument("--title", default="", help="그래프 제목 접두사 (예: 단백질 이름)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    made = 0
    for path in args.xvg:
        if not os.path.isfile(path):
            print(f"[warn] 파일 없음: {path}")
            continue
        if plot_one(path, args.outdir, args.title):
            made += 1

    if made == 0:
        sys.exit("[error] 생성된 PNG 가 없습니다.")
    print(f"[done] {made}개 PNG 생성 완료 → {args.outdir}")


if __name__ == "__main__":
    main()
