#!/usr/bin/env python3
"""Extract generator-level lepton pt, lepton eta and m_ll distributions from the showered
DY simulation.

Streams the 654 parquet skims (reading only the gen columns) and applies
  15 < m_ll < 120 GeV and lepton pt cuts (default: leading > 25 GeV via --ptLeadMin,
  subleading > 15 GeV via --ptMin; optional upper cut via --ptMax),
then plots pt, eta and m_ll for two selections:
  noetacut   no eta cut
  eta2p4     additionally |eta| < 2.4 for both leptons
Plots are normalised to expected Run 2 events (xsec_pb * lumi_pb / n_total).

  usage (inside singularity): python liv_pheno_projections.py [--jobs N] [-o/--outdir DIR]
                                  [--ptMin X] [--ptLeadMin Z] [--ptMax Y] [--linear]
"""

import argparse
import glob
import os
from concurrent.futures import ProcessPoolExecutor
from functools import partial

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from wums import output_tools, plot_tools

BASEDIR = "/ceph/submit/data/user/b/bsu17/dy_run2_showered"
SKIMDIR = os.path.join(BASEDIR, "skims")
COUNTS = os.path.join(BASEDIR, "counts.npz")

MLL_MIN, MLL_MAX = 15.0, 120.0
PT_MAX = 100.0  # plotting range only, no cut
ETA_MAX = 2.4

# name -> (eta cut or None, eta binning)
SELECTIONS = {
    "noetacut": (None, np.linspace(-10, 10, 201)),
    "eta2p4": (ETA_MAX, np.linspace(-ETA_MAX, ETA_MAX, 49)),
}
PT_BINS = np.arange(0.0, PT_MAX + 1, 1.0)
MLL_BINS = np.arange(MLL_MIN, MLL_MAX + 1, 1.0)

COLUMNS = ["gen_pt1", "gen_pt2", "gen_eta1", "gen_eta2", "gen_mass"]
VARS = ["pt", "pt1", "pt2", "eta", "eta1", "eta2", "mll"]


def bins_for(var, sel, mll_bins=MLL_BINS):
    if var.startswith("pt"):
        return PT_BINS
    if var.startswith("eta"):
        return SELECTIONS[sel][1]
    return mll_bins


def pt_cut_label(pt_min, pt_max, pt_lead_min=None):
    if pt_min is None and pt_max is None and pt_lead_min is None:
        return "none"
    parts = []
    if pt_min is not None or pt_max is not None:
        lo = "" if pt_min is None else f"{pt_min} < "
        hi = "" if pt_max is None else f" < {pt_max}"
        parts.append(f"{lo}pt{hi} GeV (both leptons)")
    if pt_lead_min is not None:
        parts.append(f"leading pt > {pt_lead_min} GeV")
    return ", ".join(parts)


def process(path, pt_min=None, pt_max=None, pt_lead_min=None, mll_bins=MLL_BINS):
    t = pq.read_table(path, columns=COLUMNS)
    t = {
        c: np.asarray(t.column(c).to_numpy(zero_copy_only=False), float)
        for c in COLUMNS
    }
    pt1, pt2, eta1, eta2, mll = (t[c] for c in COLUMNS)

    base = (mll > MLL_MIN) & (mll < MLL_MAX)
    if pt_min is not None:
        base = base & (pt1 > pt_min) & (pt2 > pt_min)
    if pt_max is not None:
        base = base & (pt1 < pt_max) & (pt2 < pt_max)
    if pt_lead_min is not None:
        base = base & (np.maximum(pt1, pt2) > pt_lead_min)
    out = {}
    for sel, (eta_cut, eta_bins) in SELECTIONS.items():
        mask = base
        if eta_cut is not None:
            mask = mask & (np.abs(eta1) < eta_cut) & (np.abs(eta2) < eta_cut)
        h = {
            "pt1": np.histogram(pt1[mask], PT_BINS)[0],
            "pt2": np.histogram(pt2[mask], PT_BINS)[0],
            "eta1": np.histogram(eta1[mask], eta_bins)[0],
            "eta2": np.histogram(eta2[mask], eta_bins)[0],
            "mll": np.histogram(mll[mask], mll_bins)[0],
        }
        h["pt"] = h["pt1"] + h["pt2"]
        h["eta"] = h["eta1"] + h["eta2"]
        h["n_pass"] = np.int64(mask.sum())
        out[sel] = h
    return out, len(mll)


def plot(H, norm, outdir, meta_info, args, mll_bins=MLL_BINS):
    import matplotlib.pyplot as plt

    for var, xlabel, logy in [
        ("pt", r"lepton $p_T$ [GeV]", True),
        ("eta", r"lepton $\eta$", False),
        ("mll", r"$m_{\ell\ell}$ [GeV]", True),
    ]:
        for sel, (eta_cut, _) in SELECTIONS.items():
            fig, ax = plt.subplots()
            rows = []
            # with custom mll binning, plot the mll spectrum per GeV
            per_gev = var == "mll" and args.mllBins is not None
            for v in (f"{var}1", f"{var}2") if var != "mll" else ("mll",):
                edges = bins_for(v, sel, mll_bins)
                y = H[sel][v] * norm
                if per_gev:
                    y = y / np.diff(edges)
                ax.stairs(y, edges, label=v)
                n_mc = int(H[sel][v].sum())
                rows.append(
                    {
                        "Histogram": v,
                        "MC events": n_mc,
                        "Yield": n_mc * norm,
                        "Uncertainty": np.sqrt(n_mc) * norm,
                    }
                )
            ax.set_xlabel(xlabel)
            ax.set_ylabel("events / GeV (Run 2)" if per_gev else "events / bin (Run 2)")
            ax.set_title(f"gen, {sel}")
            if logy and not args.linear:
                ax.set_yscale("log")
            ax.legend()

            outfile = f"{var}_{sel}"
            if args.postfix:
                outfile += f"_{args.postfix}"
            plot_tools.save_pdf_and_png(outdir, outfile, fig)
            plt.close(fig)
            output_tools.write_index_and_log(
                outdir,
                outfile,
                analysis_meta_info={
                    "Plotted histograms": pd.DataFrame(rows),
                    "Selection": {
                        "name": sel,
                        "level": "generator",
                        "mll": f"{MLL_MIN} < mll < {MLL_MAX} GeV",
                        "pt": pt_cut_label(args.ptMin, args.ptMax, args.ptLeadMin),
                        "eta": (
                            "none"
                            if eta_cut is None
                            else f"|eta| < {eta_cut} (both leptons)"
                        ),
                        "n_pass": int(H[sel]["n_pass"]),
                        "expected Run 2 events": float(H[sel]["n_pass"] * norm),
                    },
                    "Binning": {
                        "variable": var,
                        "edges": str(bins_for(var, sel, mll_bins).round(6).tolist()),
                        "normalised by bin width": per_gev,
                    },
                    "InputSample": meta_info,
                },
                args=args,
            )
    print(f"plots written to {outdir}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--glob", default="skim_*.parquet")
    ap.add_argument(
        "-o", "--outdir", default="/ceph/submit/data/user/j/jbenke/liv_pheno"
    )
    ap.add_argument(
        "-p", "--postfix", type=str, default=None, help="Postfix for output file name"
    )
    ap.add_argument(
        "--ptMin",
        type=float,
        default=15.0,
        help="min pt for both leptons [GeV] (default: 15)",
    )
    ap.add_argument(
        "--ptLeadMin",
        type=float,
        default=25.0,
        help="min pt for the leading lepton [GeV] (default: 25)",
    )
    ap.add_argument(
        "--ptMax",
        type=float,
        default=None,
        help="max pt for both leptons [GeV] (default: no cut)",
    )
    ap.add_argument(
        "--linear",
        action="store_true",
        help="linear y axis for all plots (default: log for pt, mll)",
    )
    ap.add_argument(
        "--mllBins",
        type=float,
        nargs="+",
        default=None,
        help="mll bin edges [GeV]; the mll plot is then normalised to events / GeV",
    )
    args = ap.parse_args()

    mll_bins = MLL_BINS
    if args.mllBins is not None:
        mll_bins = np.array(sorted(args.mllBins))
        if len(mll_bins) < 2:
            raise SystemExit("--mllBins needs at least two edges")
        if mll_bins[0] < MLL_MIN or mll_bins[-1] > MLL_MAX:
            print(
                f"WARNING: --mllBins extends outside the {MLL_MIN} < mll < {MLL_MAX} GeV selection, "
                "those bins will be (partially) empty"
            )

    files = sorted(glob.glob(os.path.join(SKIMDIR, args.glob)))
    if not files:
        raise SystemExit(f"no skims matched {args.glob} in {SKIMDIR}")
    print(f"processing {len(files)} skims with {args.jobs} workers", flush=True)

    H = {
        sel: {v: np.zeros(len(bins_for(v, sel, mll_bins)) - 1, np.int64) for v in VARS}
        | {"n_pass": 0}
        for sel in SELECTIONS
    }
    n_total = 0
    with ProcessPoolExecutor(args.jobs) as ex:
        work = partial(
            process,
            pt_min=args.ptMin,
            pt_max=args.ptMax,
            pt_lead_min=args.ptLeadMin,
            mll_bins=mll_bins,
        )
        for i, (out, n) in enumerate(ex.map(work, files), 1):
            for sel, h in out.items():
                for k, v in h.items():
                    H[sel][k] += v
            n_total += n
            if i % 50 == 0 or i == len(files):
                print(f"  {i}/{len(files)}", flush=True)

    counts = np.load(COUNTS)
    xsec_pb, lumi_pb = float(counts["xsec_pb"]), float(counts["lumi_pb"])
    norm = xsec_pb * lumi_pb / n_total

    for sel in SELECTIONS:
        print(
            f"{sel}: n_pass={H[sel]['n_pass']}, "
            f"expected Run 2 events={H[sel]['n_pass'] * norm:.4g}"
        )

    meta_info = {
        "base directory": BASEDIR,
        "skim directory": SKIMDIR,
        "skim files": f"{len(files)} matching {args.glob}",
        "counts file (xsec, lumi)": COUNTS,
        "n_total": n_total,
        "xsec_pb": xsec_pb,
        "lumi_pb": lumi_pb,
        "norm (events per MC event)": norm,
    }
    outdir = output_tools.make_plot_dir(args.outdir)
    plot(H, norm, outdir, meta_info, args, mll_bins)


if __name__ == "__main__":
    main()
