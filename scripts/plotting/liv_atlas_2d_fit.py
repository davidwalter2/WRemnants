#!/usr/bin/env python3
"""Expected 1D bounds on the sidereal d=4 SME coefficients (the Wilson coefficients of
Table 4 of Lunghi & Sherrill, arXiv:2412.14305) from a 2D (m_ll, y_ll) fit with rabbit.

Selection (generator level, showered DY skims of liv_pheno_projections.py), as the ATLAS
13 TeV measurement arXiv:1912.02844:
    66 <= m_ll <= 116 GeV, |eta| < 2.5 and pT >= 27 GeV for both leptons, |y_ll| < 2.5
Binning: 10 m_ll bins (5 GeV), 2 signed y_ll bins [-2.5, 0, 2.5], 4 sidereal time bins.
The sidereal time axis is made by splitting the events for the chosen luminosity (--lumi,
the first xsec * lumi generated events, so each MC event has weight 1 as for real data)
equally over the time bins. The pseudo-data are the expectation (asimov), with Poisson
uncertainties only, so the fitted coefficients are 0 and their uncertainties are the
expected bounds.

The fit model is LIV_rapidity_pheno (rabbit/param_models/sme_double_differential.py): the
relative LIV shift of each (time, m_ll, y_ll) bin from the parton-level SME basis of that
(m_ll, y_ll) bin, applied to the template yields (the lepton cuts enter through the
template, not the shift). The basis for the 10 x 2 bins is computed here and cached.

Each coefficient is fitted on its own (1D bounds, the others at zero). '--fit all' runs the
Table 2 set in its order: c_u, d_u, c_d, c_c, d_c, c_s, c_b, d_b, each for XY, XZ, YZ and
XX-YY.

  usage (inside singularity, after source setup.sh):
    python scripts/plotting/liv_atlas_2d_fit.py -o OUTDIR --fit all
    python scripts/plotting/liv_atlas_2d_fit.py -o OUTDIR --lumi 140 --fit cxyu cxzu
"""

import argparse
import glob
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from functools import partial

import hist
import numpy as np
import pyarrow.parquet as pq
from liv_pheno_fit_rapidity import int_bins_of, process
from liv_pheno_projections import COUNTS, SKIMDIR
from liv_zpole_atlas_xsec_fit import SCALE, label

from rabbit import tensorwriter
from rabbit.param_models.liv_rapidity.basis import load_or_compute_basis

MODEL = "LIV_rapidity_pheno"
# coefficients of Table 2 in its order (mode, quark), each for the directions of Table 4
TABLE2 = [
    ("c", "u"),
    ("d", "u"),
    ("c", "d"),
    ("c", "c"),
    ("d", "c"),
    ("c", "s"),
    ("c", "b"),
    ("d", "b"),
]
DIRECTIONS = ["xy", "xz", "yz", "xx"]
# all directions of the model (tx, ty: the time-space coefficients c^{TX}, c^{TY})
ALL_DIRECTIONS = DIRECTIONS + ["tx", "ty"]


def table2_coefficients(directions=DIRECTIONS):
    return [f"{m}{d}{q}" for m, q in TABLE2 for d in directions]


CHANNEL = "ch_2d"


def select_jobs(files, n_target):
    """(file, rows) jobs reading the first n_target events of the skims"""
    jobs, n_sel = [], 0
    for f in files:
        if n_sel >= n_target:
            break
        n_rows = pq.ParquetFile(f).metadata.num_rows
        take = min(n_rows, n_target - n_sel)
        jobs.append((f, take if take < n_rows else None))
        n_sel += take
    if n_sel < n_target:
        raise SystemExit(
            f"need {n_target} events, but the {len(files)} skims only have {n_sel}"
        )
    return jobs


def output_name(args):
    outname = args.outname or (
        f"liv_2d_{args.experiment}_mll{args.mllMin:g}to{args.mllMax:g}_{args.nMllBins}x{len(args.yllEdges) - 1}"
        f"_{args.nTimeBins}timeBins_{args.lumi:g}fb".replace(".", "p")
    )
    if args.postfix:
        outname += f"_{args.postfix}"
    return outname


def write_input(args):
    if args.reuseInput:
        outname = output_name(args)
        path = os.path.join(args.outdir, outname + ".hdf5")
        if os.path.exists(path):
            print(f"reusing {path}")
            return path, outname
    mll_edges = np.linspace(args.mllMin, args.mllMax, args.nMllBins + 1)
    yll_edges = np.array(args.yllEdges, dtype=float)
    axis_time = hist.axis.Regular(
        args.nTimeBins, 0, 24, name="time", overflow=False, underflow=False
    )
    axis_mll = hist.axis.Variable(
        mll_edges, name="mll", overflow=False, underflow=False
    )
    axis_yll = hist.axis.Variable(
        yll_edges, name="yll", overflow=False, underflow=False
    )

    files = sorted(glob.glob(os.path.join(SKIMDIR, args.glob)))
    if not files:
        raise SystemExit(f"no skims matched {args.glob} in {SKIMDIR}")
    counts = np.load(COUNTS)
    xsec_pb = float(counts["xsec_pb"])
    n_target = int(round(xsec_pb * args.lumi * 1e3))
    jobs = select_jobs(files, n_target)
    print(
        f"{args.lumi:g} fb^-1: reading {n_target} events from {len(jobs)} of {len(files)} skims",
        flush=True,
    )

    counts2d = np.zeros((axis_mll.size, axis_yll.size))
    n_total = 0
    work = partial(
        process,
        pt_min=args.ptMin,
        pt_lead_min=args.ptMin,
        eta_max=args.etaMax,
        mll_edges=mll_edges,
        yll_edges=yll_edges,  # events outside the m_ll and y_ll edges are dropped
    )
    with ProcessPoolExecutor(args.jobs) as ex:
        for i, (h, n) in enumerate(ex.map(work, jobs), 1):
            counts2d += h
            n_total += n
            if i % 50 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)}", flush=True)
    norm = (
        xsec_pb * args.lumi * 1e3 / n_total
    )  # = 1 up to rounding (one MC event per data event)
    n_sel = counts2d.sum() * norm
    print(f"selected {n_sel:.6g} events: {n_sel / args.nTimeBins:.6g} per time bin")
    with np.printoptions(precision=4, linewidth=150):
        print(
            "events per (m_ll, y_ll) bin per time bin:\n",
            counts2d * norm / args.nTimeBins,
        )

    # events split equally over the time bins
    h = hist.Hist(axis_time, axis_mll, axis_yll, storage=hist.storage.Weight())
    h.values()[...] = counts2d[None] * norm / args.nTimeBins
    h.variances()[...] = counts2d[None] * norm**2 / args.nTimeBins
    writer = tensorwriter.TensorWriter()
    writer.add_channel(h.axes, CHANNEL)
    writer.add_data(h.values(), CHANNEL, variances=h.values())
    writer.add_process(h, "Zll", CHANNEL, signal=True)

    outname = output_name(args)
    writer.write(
        outfolder=args.outdir,
        outfilename=outname,
        meta_data_dict={"liv_experiment": args.experiment},
    )
    path = os.path.join(args.outdir, outname + ".hdf5")
    print(f"wrote {path}")

    if not args.noBasis:
        # (m_ll, y_ll) integration bins in the order the model builds them
        int_bins = int_bins_of({a.name: a for a in h.axes})
        print(f"SME response basis for {len(int_bins)} (m_ll, y_ll) bins", flush=True)
        load_or_compute_basis(int_bins, jobs=args.jobs)
    return path, outname


def asimov_sigma(path, coeffs):
    """analytic Poisson asimov uncertainty sigma_c = (sum_i N_i s_i^2)^-1/2, s_i the relative
    shift of bin i for a unit coefficient: cross-check of the rabbit fit"""
    from rabbit import inputdata
    from rabbit.param_models import helpers as ph

    indata = inputdata.FitInputData(path)
    model = ph.load_model(MODEL, indata, *coeffs)
    names = [p.decode() for p in model.params]
    shift = model.sme.numpy() / model.sm_sigma.numpy()[None]
    nexp = np.asarray(indata.data_obs)
    return {c: 1.0 / np.sqrt(np.sum(nexp * shift[names.index(c)] ** 2)) for c in coeffs}


def run_fits(args, path, outname):
    from rabbit import io_tools

    coeffs = table2_coefficients(args.directions) if args.fit == ["all"] else args.fit
    fitdir = os.path.join(args.outdir, f"fits_{outname}")
    os.makedirs(fitdir, exist_ok=True)
    results = {}
    for c in coeffs:
        fname = f"fitresults_{c}.hdf5"
        cmd = [
            "rabbit_fit.py",
            path,
            "--paramModel",
            MODEL,
            c,
            "--noBinByBinStat",
            "-o",
            fitdir,
            "--outname",
            fname,
        ] + args.fitArgs
        print(" ".join(cmd), flush=True)
        out = subprocess.run(cmd, capture_output=True, text=True)
        if out.returncode != 0:
            print(out.stdout[-3000:], out.stderr[-3000:])
            raise SystemExit(f"fit for {c} failed")
        parms = io_tools.get_fitresult(os.path.join(fitdir, fname))["parms"].get()
        i = list(parms.axes[0]).index(c)
        results[c] = (parms.values()[i] * SCALE, np.sqrt(parms.variances()[i]) * SCALE)

    check = asimov_sigma(path, coeffs)
    print(
        f"\nexpected bounds, 2D (m_ll, y_ll) fit: {args.mllMin:g} <= m_ll <= {args.mllMax:g} GeV in "
        f"{args.nMllBins} bins, y_ll edges {args.yllEdges}, pT >= {args.ptMin:g} GeV, |eta| < {args.etaMax:g}, "
        f"{args.nTimeBins} time bins, {args.lumi:g} fb^-1, {args.experiment}, stat only"
    )
    print(
        "  1 sigma and 95% CL (1.96 sigma) of the asimov fit; 'chi2 1sig': analytic check; "
        "'95%+med': (0.674 + 1.96) sigma, the median 95% bound of pseudo-experiments"
    )
    print(
        f"{'coefficient':13s} {'1 sigma':>10s} {'95% CL':>10s} {'95%+med':>10s} {'chi2 1sig':>10s}"
    )
    for c, (v, e) in results.items():
        print(
            f"{label(c):13s} {e:10.2e} {1.96 * e:10.2e} {(0.6745 + 1.96) * e:10.2e} {check[c]:10.2e}"
        )
    # merge with the bounds of earlier runs into the same fit directory (new fits replace old ones)
    bounds_file = os.path.join(fitdir, "bounds.npz")
    merged = {}
    if os.path.exists(bounds_file):
        old = np.load(bounds_file)
        for c, v, e, a in zip(
            old["coefficients"], old["values"], old["sigma"], old["sigma_analytic"]
        ):
            merged[str(c)] = (v, e, a)
    for c, (v, e) in results.items():
        merged[c] = (v, e, check[c])
    np.savez(
        bounds_file,
        coefficients=np.array(list(merged)),
        values=np.array([m[0] for m in merged.values()]),
        sigma=np.array([m[1] for m in merged.values()]),
        sigma_analytic=np.array([m[2] for m in merged.values()]),
    )
    print(f"saved {os.path.join(fitdir, 'bounds.npz')}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("-o", "--outdir", default=os.getcwd())
    ap.add_argument("--outname", default=None)
    ap.add_argument("-p", "--postfix", default=None)
    ap.add_argument(
        "--jobs", type=int, default=16, help="workers for the skims and the SME basis"
    )
    ap.add_argument("--glob", default="skim_*.parquet")
    ap.add_argument(
        "--lumi",
        type=float,
        default=36.1,
        help="integrated luminosity [fb^-1] (default: ATLAS 2015-16)",
    )
    ap.add_argument(
        "--experiment", default="ATLAS", type=str.upper, choices=["CMS", "ATLAS"]
    )
    ap.add_argument("--nTimeBins", type=int, default=4)
    ap.add_argument("--mllMin", type=float, default=66.0)
    ap.add_argument("--mllMax", type=float, default=116.0)
    ap.add_argument("--nMllBins", type=int, default=10)
    ap.add_argument(
        "--yllEdges",
        type=float,
        nargs="+",
        default=[-2.5, 0.0, 2.5],
        help="signed y_ll bin edges",
    )
    ap.add_argument(
        "--ptMin", type=float, default=27.0, help="min pT of both leptons [GeV]"
    )
    ap.add_argument(
        "--etaMax", type=float, default=2.5, help="max |eta| of both leptons"
    )
    ap.add_argument(
        "--noBasis",
        action="store_true",
        help="do not precompute the SME response basis",
    )
    ap.add_argument(
        "--fit",
        nargs="*",
        default=None,
        help="1D rabbit fit per coefficient ('all' for the Table 2 set) and print the bounds",
    )
    ap.add_argument(
        "--directions",
        nargs="+",
        default=DIRECTIONS,
        choices=ALL_DIRECTIONS,
        help="directions of the coefficients fitted with '--fit all'",
    )
    ap.add_argument(
        "--reuseInput",
        action="store_true",
        help="use the existing fit input if it exists",
    )
    ap.add_argument(
        "--fitArgs",
        nargs=argparse.REMAINDER,
        default=[],
        help="extra rabbit_fit.py arguments",
    )
    args = ap.parse_args()

    path, outname = write_input(args)
    if args.fit is not None:
        run_fits(args, path, outname)
    else:
        print(
            f"fit with: rabbit_fit.py {path} --paramModel {MODEL} cxyu --noBinByBinStat"
        )


if __name__ == "__main__":
    main()
