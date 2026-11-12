#!/usr/bin/env python3
"""Extract parton-level lepton pt, lepton eta and m_ll distributions from the MadGraph
DY LHE file (pp -> l+l- + 0/1 jet, MLM matched, before parton shower / Delphes).

Streams the LHE file(s) in BASEDIR, takes the two final-state charged leptons of each
event (lepton 1 = higher pt), applies
  15 < m_ll < 120 GeV  and  15 < pt < 80 GeV for both leptons,
then plots pt, eta and m_ll for two selections:
  noetacut   no eta cut
  eta2p4     additionally |eta| < 2.4 for both leptons
Plots are normalised to expected Run 2 events (xsec_pb * lumi_pb / n_total), with the
cross section read from the LHE <init> block (summed over all subprocesses).
Note: the MLM matching veto is applied in the shower, so LHE-level 0j+1j is not yet merged.

  usage (inside singularity): python liv_pheno_projections_lhe.py [--jobs N] [-o/--outdir DIR]
"""

import argparse
import glob
import gzip
import os
import re
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from wums import output_tools, plot_tools

# BASEDIR = "/ceph/submit/data/user/b/bsu17/dy_run2"  # old sample: events_*.lhe.gz
BASEDIR = "/work/submit/bsu17/MG5/jets/Events/run_01"
# Run 2 luminosity used for the showered sample (dy_run2_showered/counts.npz)
LUMI_PB = 140168.88850786575

MLL_MIN, MLL_MAX = 15.0, 120.0
PT_MIN, PT_MAX = 15.0, 80.0
ETA_MAX = 2.4

# name -> (eta cut or None, eta binning)
SELECTIONS = {
    "noetacut": (None, np.linspace(-10, 10, 201)),
    "eta2p4": (ETA_MAX, np.linspace(-ETA_MAX, ETA_MAX, 49)),
}
PT_BINS = np.arange(PT_MIN, PT_MAX + 1, 1.0)
MLL_BINS = np.arange(MLL_MIN, MLL_MAX + 1, 1.0)

VARS = ["pt", "pt1", "pt2", "eta", "eta1", "eta2", "mll"]

# LHE particle line: pid status mother1 mother2 color1 color2 px py pz E m ...
# keep final-state (status 1) charged leptons
LEPTON_RE = re.compile(
    rb"^\s*(-?1[135])\s+1\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", re.M
)
INIT_RE = re.compile(rb"<init>\s*\n(.*?)</init>", re.S)


def bins_for(var, sel):
    if var.startswith("pt"):
        return PT_BINS
    if var.startswith("eta"):
        return SELECTIONS[sel][1]
    return MLL_BINS


def read_leptons(path):
    with gzip.open(path, "rb") as f:
        data = f.read()
    # <init>: beam line ending in NPRUP, then one "xsec xerr xmax lprup" line per subprocess
    init = INIT_RE.search(data).group(1).strip().splitlines()
    nprup = int(init[0].split()[-1])
    xsec_pb = sum(float(l.split()[0]) for l in init[1 : 1 + nprup])
    n_events = data.count(b"<event>")

    lep = np.array(LEPTON_RE.findall(data), dtype=float)
    if len(lep) != 2 * n_events:
        raise RuntimeError(f"{path}: found {len(lep)} leptons for {n_events} events")
    lep = lep.reshape(n_events, 2, 5)  # (event, lepton, [pid, px, py, pz, E])
    return lep, n_events, xsec_pb


def process(path):
    lep, n_events, xsec_pb = read_leptons(path)
    px, py, pz, E = (lep[..., i] for i in range(1, 5))

    pt = np.hypot(px, py)
    eta = np.arcsinh(pz / pt)
    mll = np.sqrt(
        np.maximum(E.sum(1) ** 2 - px.sum(1) ** 2 - py.sum(1) ** 2 - pz.sum(1) ** 2, 0)
    )

    # order leptons by pt
    order = np.argsort(-pt, axis=1)
    pt = np.take_along_axis(pt, order, 1)
    eta = np.take_along_axis(eta, order, 1)
    pt1, pt2, eta1, eta2 = pt[:, 0], pt[:, 1], eta[:, 0], eta[:, 1]

    base = (
        (mll > MLL_MIN)
        & (mll < MLL_MAX)
        & (pt1 > PT_MIN)
        & (pt1 < PT_MAX)
        & (pt2 > PT_MIN)
        & (pt2 < PT_MAX)
    )
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
            "mll": np.histogram(mll[mask], MLL_BINS)[0],
        }
        h["pt"] = h["pt1"] + h["pt2"]
        h["eta"] = h["eta1"] + h["eta2"]
        h["n_pass"] = np.int64(mask.sum())
        out[sel] = h
    return out, n_events, xsec_pb


def plot(H, norm, outdir, meta_info, args):
    import matplotlib.pyplot as plt

    for var, xlabel, logy in [
        ("pt", r"lepton $p_T$ [GeV]", True),
        ("eta", r"lepton $\eta$", False),
        ("mll", r"$m_{\ell\ell}$ [GeV]", True),
    ]:
        for sel, (eta_cut, _) in SELECTIONS.items():
            fig, ax = plt.subplots()
            rows = []
            for v in (f"{var}1", f"{var}2") if var != "mll" else ("mll",):
                ax.stairs(H[sel][v] * norm, bins_for(v, sel), label=v)
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
            ax.set_ylabel("events / bin (Run 2)")
            ax.set_title(f"LHE, {sel}")
            if logy:
                ax.set_yscale("log")
            ax.legend()

            outfile = f"{var}_{sel}"
            plot_tools.save_pdf_and_png(outdir, outfile, fig)
            plt.close(fig)
            output_tools.write_index_and_log(
                outdir,
                outfile,
                analysis_meta_info={
                    "Plotted histograms": pd.DataFrame(rows),
                    "Selection": {
                        "name": sel,
                        "level": "LHE (parton level, before shower)",
                        "mll": f"{MLL_MIN} < mll < {MLL_MAX} GeV",
                        "pt": f"{PT_MIN} < pt < {PT_MAX} GeV (both leptons)",
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
                        "edges": str(bins_for(var, sel).round(6).tolist()),
                    },
                    "InputSample": meta_info,
                },
                args=args,
            )
    print(f"plots written to {outdir}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--glob", default="unweighted_events.lhe.gz")
    ap.add_argument(
        "--lumi", type=float, default=LUMI_PB, help="integrated luminosity in pb^-1"
    )
    ap.add_argument(
        "-o", "--outdir", default="/ceph/submit/data/user/j/jbenke/liv_pheno_lhe"
    )
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(BASEDIR, args.glob)))
    if not files:
        raise SystemExit(f"no LHE files matched {args.glob} in {BASEDIR}")
    print(f"processing {len(files)} LHE files with {args.jobs} workers", flush=True)

    H = {
        sel: {v: np.zeros(len(bins_for(v, sel)) - 1, np.int64) for v in VARS}
        | {"n_pass": 0}
        for sel in SELECTIONS
    }
    n_total = 0
    xsecs = []
    with ProcessPoolExecutor(args.jobs) as ex:
        for i, (out, n, xsec) in enumerate(ex.map(process, files), 1):
            for sel, h in out.items():
                for k, v in h.items():
                    H[sel][k] += v
            n_total += n
            xsecs.append(xsec)
            if i % 50 == 0 or i == len(files):
                print(f"  {i}/{len(files)}", flush=True)

    # unweighted events: each file's sample is equally sized, so average the xsec
    xsec_pb = float(np.mean(xsecs))
    norm = xsec_pb * args.lumi / n_total

    for sel in SELECTIONS:
        print(
            f"{sel}: n_pass={H[sel]['n_pass']}, "
            f"expected Run 2 events={H[sel]['n_pass'] * norm:.4g}"
        )

    meta_info = {
        "base directory": BASEDIR,
        "LHE files": f"{len(files)} matching {args.glob}",
        "generator": "MadGraph5_aMC@NLO 2.9.27, pp -> l+l- (+j), MLM xqcut=30, ptl > 8 GeV, mmll > 12 GeV, 13 TeV",
        "n_total": n_total,
        "xsec_pb (from LHE <init>)": xsec_pb,
        "lumi_pb": args.lumi,
        "norm (events per MC event)": norm,
    }
    outdir = output_tools.make_plot_dir(args.outdir)
    plot(H, norm, outdir, meta_info, args)


if __name__ == "__main__":
    main()
