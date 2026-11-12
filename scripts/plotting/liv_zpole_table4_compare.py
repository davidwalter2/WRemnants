#!/usr/bin/env python3
"""Compare the expected 1D Z-pole bounds of liv_zpole_atlas_xsec_fit.py with Table 4 of
Lunghi & Sherrill, arXiv:2412.14305.

Reads bounds.npz of the fits with --uncorrScaling none and th, lumi, syst all 100% correlated
between the time bins (the last column of Table 4), for dsigma/dm_ll at m_Z integrated over
all y_ll (the paper's observable) and over |y_ll| < 2.5. Our bounds are the asimov 95% CL,
1.96 sigma. Prints the table and saves it as text and csv.

  usage: python scripts/plotting/liv_zpole_table4_compare.py -i liv_zpole_table4
"""

import argparse
import csv
import os

import numpy as np

# Table 4 of arXiv:2412.14305 (transcribed from the arXiv HTML, v1): columns delta_th,
# (delta_th, delta_lumi), (delta_th, delta_lumi, delta_syst) taken 100% correlated; only the
# last one is compared
PAPER = {
    "cxyu": (2.8e-4, 1.0e-4, 3.1e-6),
    "cxzu": (9.4e-4, 3.4e-4, 1.1e-5),
    "cyzu": (9.3e-4, 3.3e-4, 1.0e-5),
    "cxxu": (2.0e-3, 7.2e-4, 2.3e-5),
    "cxyd": (1.2e-2, 4.4e-3, 1.4e-4),
    "cxzd": (4.2e-2, 1.5e-2, 4.7e-4),
    "cyzd": (4.1e-2, 1.5e-2, 4.6e-4),
    "cxxd": (8.9e-2, 3.2e-2, 1.0e-3),
    "dxyu": (3.6e-4, 1.3e-4, 4.0e-6),
    "dxzu": (1.2e-3, 4.3e-4, 1.4e-5),
    "dyzu": (1.2e-3, 4.3e-4, 1.3e-5),
    "dxxu": (2.6e-3, 9.3e-4, 2.9e-5),
    "cxyc": (5.3e-3, 1.9e-3, 6.0e-5),
    "cxzc": (1.8e-2, 6.5e-3, 2.0e-4),
    "cyzc": (1.8e-2, 6.4e-3, 2.0e-4),
    "cxxc": (3.8e-2, 1.4e-2, 4.4e-4),
    "cxys": (1.1e-1, 3.9e-2, 1.2e-3),
    "cxzs": (3.7e-1, 1.3e-1, 4.2e-3),
    "cyzs": (3.6e-1, 1.3e-1, 4.0e-3),
    "cxxs": (7.9e-1, 2.8e-1, 8.9e-3),
    "dxyc": (5.9e-3, 2.1e-3, 6.7e-5),
    "dxzc": (2.0e-2, 7.2e-3, 2.3e-4),
    "dyzc": (2.0e-2, 7.0e-3, 2.2e-4),
    "dxxc": (4.3e-2, 1.5e-2, 4.8e-4),
    "cxyb": (3.6e-2, 1.3e-2, 4.1e-4),
    "cxzb": (1.2e-1, 4.4e-2, 1.4e-3),
    "cyzb": (1.2e-1, 4.3e-2, 1.3e-3),
    "cxxb": (2.6e-1, 9.3e-2, 2.9e-3),
    "dxyb": (3.8e-2, 1.4e-2, 4.3e-4),
    "dxzb": (1.3e-1, 4.7e-2, 1.5e-3),
    "dyzb": (1.3e-1, 4.6e-2, 1.4e-3),
    "dxxb": (2.8e-1, 9.9e-2, 3.1e-3),
}
CORR = "th_lumi_syst"
YLL = [("all", "all y_ll"), ("2p5", "|y_ll|<2.5")]


def label(c):
    d = "XX-YY" if c[1:3] == "xx" else c[1:3].upper()
    return f"{c[0]}_{c[3]}^{d}"


def load(indir, yll, corr, ntime):
    name = f"liv_zpole_atlas_mZ_yll{yll}_{ntime}timeBins_corr_{corr}_uncorrScalingNone"
    f = np.load(os.path.join(indir, f"fits_{name}", "bounds.npz"))
    return {
        str(c): (s, s2)
        for c, s, s2 in zip(f["coefficients"], f["sigma"], f["sigma_chi2"])
    }


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("-i", "--indir", default="liv_zpole_table4")
    ap.add_argument("--nTimeBins", type=int, default=4)
    ap.add_argument(
        "--bounds2d",
        default="liv_atlas_2d/fits_liv_2d_ATLAS_mll66to116_10x2_4timeBins_36p1fb/bounds.npz",
        help="bounds.npz of liv_atlas_2d_fit.py (2D m_ll, y_ll fiducial fit, stat only), compared with "
        "the all-y_ll 1D bounds; '' to skip",
    )
    ap.add_argument(
        "-o",
        "--outname",
        default="table4_comparison",
        help="written to INDIR/OUTNAME.{txt,csv}",
    )
    args = ap.parse_args()

    res = {y: load(args.indir, y, CORR, args.nTimeBins) for y, _ in YLL}
    coeffs = list(PAPER)
    res2d = {}
    if args.bounds2d:
        f = np.load(args.bounds2d)
        res2d = {str(c): s for c, s in zip(f["coefficients"], f["sigma"])}

    rows, lines = [], []
    lines.append(
        "Expected 1D bounds, dsigma/dm_ll at m_Z, ATLAS 36.1 fb^-1, 4 sidereal time bins"
    )
    lines.append(
        "th, lumi, syst 100% correlated between time bins; uncorrelated (stat) uncertainty not "
        "scaled by sqrt(n) (--uncorrScaling none)"
    )
    lines.append(
        "ours: asimov 95% CL = 1.96 sigma. paper: Table 4 of arXiv:2412.14305, last column"
    )
    if res2d:
        lines.append(
            f"2D: {args.bounds2d}: 66-116 GeV x 10 m_ll, 2 signed y_ll bins, pT > 27, |eta|, |y_ll| < 2.5, "
            "stat only (Poisson, events split over the time bins), asimov 1.96 sigma"
        )
    lines.append("")
    head2d = f" {'2D':>9s} {'2D/all':>7s}" if res2d else ""
    lines.append(
        f"{'coefficient':12s} {'paper':>9s} {'all y_ll':>9s} {'paper/all':>9s} "
        f"{'|y|<2.5':>9s} {'2.5/all':>8s}" + head2d
    )
    for c in coeffs:
        p = PAPER[c][2]
        s_all, chi_all = res["all"][c]
        s_25, chi_25 = res["2p5"][c]
        b_all, b_25 = 1.96 * s_all, 1.96 * s_25
        col2d = ""
        extra = {}
        if res2d:
            b_2d = 1.96 * res2d[c]
            col2d = f" {b_2d:9.2e} {b_2d / b_all:7.1f}"
            extra = dict(
                sigma_2d=res2d[c], bound95_2d=b_2d, ratio_2d_over_yllall=b_2d / b_all
            )
        lines.append(
            f"{label(c):12s} {p:9.2e} {b_all:9.2e} {p / b_all:9.2f} {b_25:9.2e} {b_25 / b_all:8.1f}"
            + col2d
        )
        rows.append(
            dict(
                coefficient=c,
                paper_95=p,
                sigma_yllall=s_all,
                bound95_yllall=b_all,
                sigma_chi2_yllall=chi_all,
                sigma_yll2p5=s_25,
                bound95_yll2p5=b_25,
                sigma_chi2_yll2p5=chi_25,
                ratio_paper_over_yllall=p / b_all,
                ratio_yll2p5_over_yllall=b_25 / b_all,
                **extra,
            )
        )
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(args.indir, args.outname + ".txt"), "w") as f:
        f.write(text + "\n")
    with open(os.path.join(args.indir, args.outname + ".csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nsaved {os.path.join(args.indir, args.outname)}.txt and .csv")


if __name__ == "__main__":
    main()
