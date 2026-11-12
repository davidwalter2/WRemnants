#!/usr/bin/env python3
"""Rabbit input for a generator level 2D (m_ll, y_ll) LIV fit to the showered DY simulation
(the skims of liv_pheno_projections.py), fit with the double differential SME model
(--paramModel LIV_rapidity_pheno, rabbit/param_models/sme_double_differential.py).

Selection: 15 < m_ll < 120 GeV, leading lepton pt > 25 GeV, subleading pt > 15 GeV,
|y_ll| < 2.4, optional |eta| < X for both leptons (--etaMax, default: no cut beyond the skims).
Binning: 1 GeV in m_ll, 0.1 in y_ll. The (m_ll, y_ll) histogram, normalised to expected
Run 2 events, is split equally over the sidereal time bins (24 by default, --nTimeBins), giving one channel
(time, mll, yll). The data are set to the expectation (asimov).

The SME response basis for the 105 x 48 (m_ll, y_ll) bins is precomputed here in parallel
and cached in rabbit/param_models/precomputed_sigma, where the model picks it up
(computing it serially in the fit would take hours).

With --singleDiff the y_ll axis and cut are dropped: one (time, mll) channel, fit with
--paramModel LIV_mll_pheno, which integrates the SME model over the full rapidity range.

  usage (inside singularity, after source setup.sh):
    python liv_pheno_fit_rapidity.py [--jobs N] [-o OUTDIR]
    rabbit_fit.py OUTDIR/liv_pheno_rapidity.hdf5 --paramModel LIV_rapidity_pheno cxxu cxyu cxzu cyzu \\
        --saveHists --computeHistErrors -m Project ch_pheno time --postfix cu_pheno
    python liv_pheno_fit_rapidity.py --singleDiff [--jobs N] [-o OUTDIR]
    rabbit_fit.py OUTDIR/liv_pheno_mll.hdf5 --paramModel LIV_mll_pheno cxxu cxyu cxzu cyzu \\
        --saveHists --computeHistErrors -m Project ch_pheno time --postfix cu_pheno_mll

With --experiment ATLAS the input is identical (the time axis is read as the local
sidereal time of the experiment), but it is written as liv_pheno_*_ATLAS.hdf5 with
'liv_experiment' in the meta data, so the model uses the ATLAS beam orientation
(functions.EXPERIMENTS: beam antiparallel to CMS, latitude 46.24 deg). The orientation can
also be forced in the fit, e.g. --paramModel LIV_rapidity_pheno_ATLAS cxxu ... or
--paramModel LIV_rapidity_pheno experiment=ATLAS cxxu ...
    python liv_pheno_fit_rapidity.py --experiment ATLAS [--jobs N] [-o OUTDIR]
    rabbit_fit.py OUTDIR/liv_pheno_rapidity_ATLAS.hdf5 --paramModel LIV_rapidity_pheno cxxu cxyu cxzu cyzu \\
        --saveHists --computeHistErrors -m Project ch_pheno time --postfix cu_pheno_ATLAS

With --zpole the input is the Z-pole cross section per sidereal time bin, as in the d=4
sidereal analysis of Lunghi & Sherrill, arXiv:2412.14305 (Table 4, ATLAS 13 TeV, 36.1 fb^-1):
a single (time, mll) channel ch_z, structured like --singleDiff but with one mll bin (--zpoleMll, default mZ +- 0.5 GeV, close to the
narrow-width Z pole of the paper; the ATLAS measurement behind its uncertainties,
arXiv:1912.02844, is fiducial in 66-116 GeV, pT > 27 GeV, |eta| < 2.5) and one yll bin, |y_ll| < --zpoleYllMax ('all' for the full
rapidity range, or a number, e.g. 2.4) applied as a cut, without a yll axis. Fit with LIV_mll_pheno, which
integrates the SME template over |y| < --zpoleYllMax (read from 'liv_yll_max' in the meta data, or forced
with a 'yllMax=X' argument of --paramModel). The paper is inclusive in the lepton kinematics, use
--ptMin 0 --ptLeadMin 0 for the closest comparison, or --atlasFiducial for the lepton cuts of the
ATLAS measurement (pT > 27 GeV and |eta| < 2.5 for both leptons; the mll window stays --zpoleMll,
use --zpoleMll 66 116 for the ATLAS one). Fit one coefficient at a time for 1D limits.
    python liv_pheno_fit_rapidity.py --zpole --zpoleYllMax all --experiment ATLAS --lumi 36.1 --nTimeBins 4 \\
        --ptMin 0 --ptLeadMin 0
    rabbit_fit.py liv_pheno_zpole_yllall_noptcut_ATLAS_4timeBins_36_1fb.hdf5 --paramModel LIV_mll_pheno cxxu \\
        -m Project ch_z time -p zpole_yllall_cxxu
"""

import argparse
import glob
import os
from concurrent.futures import ProcessPoolExecutor
from functools import partial

import hist
import numpy as np
import pyarrow.parquet as pq
from liv_pheno_projections import COUNTS, SKIMDIR

from rabbit import tensorwriter
from rabbit.param_models.liv_rapidity.basis import (
    YLL_FULL,
    basis_cache_file,
    load_or_compute_basis,
)
from rabbit.param_models.liv_rapidity.functions import m_Z

MU_MASS = 0.1056583755
COLUMNS = [
    "gen_pt1",
    "gen_pt2",
    "gen_eta1",
    "gen_eta2",
    "gen_phi1",
    "gen_phi2",
    "gen_mass",
]

axis_mll = hist.axis.Regular(
    105, 15.0, 120.0, name="mll", overflow=False, underflow=False
)
axis_yll = hist.axis.Regular(48, -2.4, 2.4, name="yll", overflow=False, underflow=False)


def dilepton_rapidity(pt1, pt2, eta1, eta2, phi1, phi2):
    px = pt1 * np.cos(phi1) + pt2 * np.cos(phi2)
    py = pt1 * np.sin(phi1) + pt2 * np.sin(phi2)
    pz = pt1 * np.sinh(eta1) + pt2 * np.sinh(eta2)
    e = np.sqrt((pt1 * np.cosh(eta1)) ** 2 + MU_MASS**2) + np.sqrt(
        (pt2 * np.cosh(eta2)) ** 2 + MU_MASS**2
    )
    return 0.5 * np.log((e + pz) / (e - pz))


def process(
    job, pt_min, pt_lead_min, eta_max=None, mll_edges=None, yll_edges=None, yll_max=None
):
    path, max_rows = job
    t = pq.read_table(path, columns=COLUMNS)
    if max_rows is not None:
        t = t.slice(0, max_rows)
    t = {
        c: np.asarray(t.column(c).to_numpy(zero_copy_only=False), float)
        for c in COLUMNS
    }
    pt1, pt2, eta1, eta2, phi1, phi2, mll = (t[c] for c in COLUMNS)

    mask = (np.minimum(pt1, pt2) > pt_min) & (np.maximum(pt1, pt2) > pt_lead_min)
    if eta_max is not None:
        mask &= (np.abs(eta1) < eta_max) & (np.abs(eta2) < eta_max)
    if yll_edges is None:
        # no rapidity axis: optional |y_ll| < yll_max cut, out of range m_ll events are
        # dropped by histogram
        if yll_max is not None:
            yll = dilepton_rapidity(pt1, pt2, eta1, eta2, phi1, phi2)
            mask &= np.abs(yll) < yll_max
        h, _ = np.histogram(mll[mask], bins=mll_edges)
        return h, len(mll)
    yll = dilepton_rapidity(
        pt1[mask], pt2[mask], eta1[mask], eta2[mask], phi1[mask], phi2[mask]
    )
    # out of range events (m_ll, y_ll) are dropped by histogram2d
    h, _, _ = np.histogram2d(mll[mask], yll, bins=[mll_edges, yll_edges])
    return h, len(mll)


def int_bins_of(axes, yll_max=YLL_FULL):
    """(mll, rapidity) integration bins in the order LIV_rapidity builds them,
    without a yll axis |y| < yll_max (default: the full rapidity range of LIV_mll_pheno)
    """
    mll = axes["mll"].edges
    y = axes["yll"].edges if "yll" in axes else np.array([-yll_max, yll_max])
    return [
        (float(mlo), float(mhi), float(ylo), float(yhi))
        for mlo, mhi in zip(mll[:-1], mll[1:])
        for ylo, yhi in zip(y[:-1], y[1:])
    ]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--glob", default="skim_*.parquet")
    ap.add_argument(
        "-o",
        "--outdir",
        default=os.path.dirname(os.path.abspath(__file__)),
        help="output directory (default: the directory of this script, scripts/plotting)",
    )
    ap.add_argument(
        "--outname",
        default=None,
        help="default: liv_pheno_rapidity (liv_pheno_mll with --singleDiff)",
    )
    ap.add_argument(
        "--ptMin", type=float, default=15.0, help="min pt for both leptons [GeV]"
    )
    ap.add_argument(
        "--ptLeadMin",
        type=float,
        default=25.0,
        help="min pt for the leading lepton [GeV]",
    )
    ap.add_argument(
        "--etaMax",
        type=float,
        default=None,
        help="max |eta| for both leptons (default: no additional cut)",
    )
    ap.add_argument(
        "--singleDiff",
        action="store_true",
        help="fit (time, mll) only, without the y_ll axis and cut",
    )
    ap.add_argument(
        "-p", "--postfix", default=None, help="postfix appended to the output file name"
    )
    ap.add_argument(
        "--experiment",
        default="CMS",
        type=str.upper,
        choices=["CMS", "ATLAS"],
        help="beam orientation used by the SME model in the fit (stored in the meta data)",
    )
    ap.add_argument(
        "--nTimeBins",
        type=int,
        default=24,
        help="number of equal sidereal time bins over 24 h (4 reproduces the binning of "
        "Lunghi & Sherrill, arXiv:2412.14305, Table 4); output name gets a _<N>timeBins suffix if != 24",
    )
    ap.add_argument(
        "--lumi",
        type=float,
        default=None,
        help="integrated luminosity [fb^-1]: only the first xsec*lumi generated events of the skims are read "
        "(whole skims in --glob order, the last one truncated), so each MC event has weight 1 as for real "
        "data (default: all skims, normalised to the Run 2 lumi of the counts file); "
        "output name gets a _<L>fb suffix if set (decimal point as '_', e.g. _36_1fb)",
    )
    ap.add_argument(
        "--zpole",
        action="store_true",
        help="Z-pole input of arXiv:2412.14305: one channel ch_z with a single (mll, yll) bin per time bin",
    )
    ap.add_argument(
        "--zpoleMll",
        type=float,
        nargs=2,
        default=[m_Z - 0.5, m_Z + 0.5],
        metavar=("LO", "HI"),
        help="mll window [GeV] for --zpole (default: mZ +- 0.5, the narrow-width Z pole of arXiv:2412.14305)",
    )
    ap.add_argument(
        "--zpoleYllMax",
        default=None,
        help="|y_ll| range for --zpole, a number or 'all' (full rapidity range). Default: 'all', "
        "2.5 with --atlasFiducial (the lepton acceptance, so the SME prediction is integrated "
        "over the rapidities of the selected events)",
    )
    ap.add_argument(
        "--atlasFiducial",
        action="store_true",
        help="lepton cuts of the ATLAS 13 TeV measurement (arXiv:1912.02844): pT > 27 GeV and "
        "|eta| < 2.5 for both leptons (overrides --ptMin, --ptLeadMin, --etaMax)",
    )
    ap.add_argument(
        "--noBasis",
        action="store_true",
        help="do not precompute the SME response basis",
    )
    args = ap.parse_args()
    if args.atlasFiducial:
        args.ptMin, args.ptLeadMin, args.etaMax = 27.0, 27.0, 2.5
    if args.zpoleYllMax is None:
        args.zpoleYllMax = "2.5" if args.atlasFiducial else "all"
    if args.zpole:
        if args.singleDiff:
            raise SystemExit("--zpole and --singleDiff are exclusive")
        zpole_ymax = (
            YLL_FULL if args.zpoleYllMax.lower() == "all" else float(args.zpoleYllMax)
        )
    if args.outname is None:
        if args.zpole:
            yll_tag = (
                "all"
                if args.zpoleYllMax.lower() == "all"
                else f"{zpole_ymax:g}".replace(".", "p")
            )
            args.outname = f"liv_pheno_zpole_yll{yll_tag}"
            if args.zpoleMll != [m_Z - 0.5, m_Z + 0.5]:
                args.outname += (
                    f"_mll{args.zpoleMll[0]:g}to{args.zpoleMll[1]:g}".replace(".", "p")
                )
            if args.atlasFiducial:
                args.outname += "_atlasfid"
            elif args.ptMin <= 0 and args.ptLeadMin <= 0:
                args.outname += "_noptcut"
        else:
            args.outname = "liv_pheno_mll" if args.singleDiff else "liv_pheno_rapidity"
        if args.experiment != "CMS":
            args.outname += f"_{args.experiment}"
    if args.nTimeBins != 24:
        args.outname += f"_{args.nTimeBins}timeBins"
    if args.lumi is not None:
        args.outname += f"_{args.lumi:g}fb".replace(".", "_")
    if args.postfix:
        args.outname += f"_{args.postfix}"
    axis_time = hist.axis.Regular(
        args.nTimeBins, 0, 24, name="time", overflow=False, underflow=False
    )
    yll_max = None  # |y_ll| cut without a yll axis (None: no cut)
    if args.zpole:
        # same structure as --singleDiff: a (time, mll) channel, here with a single mll bin;
        # the model integrates |y| < zpole_ymax (written to the meta data as liv_yll_max)
        mll_edges = list(args.zpoleMll)
        yll_edges = None
        if zpole_ymax < YLL_FULL:
            yll_max = zpole_ymax
        kin_axes = [
            hist.axis.Variable(mll_edges, name="mll", overflow=False, underflow=False)
        ]
    else:
        kin_axes = [axis_mll] if args.singleDiff else [axis_mll, axis_yll]
        mll_edges = axis_mll.edges
        yll_edges = None if args.singleDiff else axis_yll.edges

    files = sorted(glob.glob(os.path.join(SKIMDIR, args.glob)))
    if not files:
        raise SystemExit(f"no skims matched {args.glob} in {SKIMDIR}")
    counts = np.load(COUNTS)
    xsec_pb = float(counts["xsec_pb"])
    if args.lumi is None:
        lumi_pb = float(counts["lumi_pb"])
        jobs = [(f, None) for f in files]
    else:
        lumi_pb = args.lumi * 1e3
        n_target = int(round(xsec_pb * lumi_pb))
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
                f"--lumi {args.lumi:g} needs {n_target} events, but the {len(files)} skims only have {n_sel} "
                f"({n_sel / xsec_pb / 1e3:.4g} fb^-1)"
            )
        print(
            f"--lumi {args.lumi:g} fb^-1: reading {n_target} events from {len(jobs)} of {len(files)} skims"
        )
    print(f"processing {len(jobs)} skims with {args.jobs} workers", flush=True)

    counts2d = np.zeros([a.size for a in kin_axes])
    n_total = 0
    with ProcessPoolExecutor(args.jobs) as ex:
        work = partial(
            process,
            pt_min=args.ptMin,
            pt_lead_min=args.ptLeadMin,
            eta_max=args.etaMax,
            mll_edges=mll_edges,
            yll_edges=yll_edges,
            yll_max=yll_max,
        )
        for i, (h, n) in enumerate(ex.map(work, jobs), 1):
            counts2d += h
            n_total += n
            if i % 50 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)}", flush=True)

    # = 1 with --lumi (one MC event per expected data event)
    norm = xsec_pb * lumi_pb / n_total
    print(
        f"n_pass={int(counts2d.sum())}, expected events for {lumi_pb / 1e3:.4g} fb^-1={counts2d.sum() * norm:.4g}"
    )

    if args.zpole:
        channels = {"ch_z": counts2d}
        n_exp = counts2d.sum() * norm
        print(
            f"  ch_z: mll {mll_edges}, |yll| < {zpole_ymax:.4g}: {n_exp:.6g} events, "
            f"stat. {1 / np.sqrt(n_exp):.3g} (total), {np.sqrt(args.nTimeBins / n_exp):.3g} (per time bin)"
        )
    else:
        channels = {"ch_pheno": counts2d}

    writer = tensorwriter.TensorWriter()
    int_bins = []
    for name, c in channels.items():
        axes = kin_axes
        # events split equally over the time bins (each gets 1/nTimeBins of the MC events)
        h = hist.Hist(axis_time, *axes, storage=hist.storage.Weight())
        h.values()[...] = c[None] * norm / args.nTimeBins
        h.variances()[...] = c[None] * norm**2 / args.nTimeBins
        writer.add_channel(h.axes, name)
        writer.add_data(h.values(), name, variances=h.values())
        writer.add_process(h, "Zmumu", name, signal=True)
        # in the order the model builds them (channel by channel)
        for b in int_bins_of(
            {a.name: a for a in h.axes}, zpole_ymax if args.zpole else YLL_FULL
        ):
            if b not in int_bins:
                int_bins.append(b)
    meta = {"liv_experiment": args.experiment}
    if args.zpole:
        meta["liv_yll_max"] = zpole_ymax
    writer.write(
        outfolder=args.outdir,
        outfilename=args.outname,
        meta_data_dict=meta,
    )
    model = "LIV_mll_pheno" if (args.singleDiff or args.zpole) else "LIV_rapidity_pheno"
    print(
        f"fit with: rabbit_fit.py {os.path.join(args.outdir, args.outname)}.hdf5 --paramModel {model} "
        f"cxxu cxyu cxzu cyzu ...  ({args.experiment} beam orientation from the meta data)"
    )

    if not args.noBasis:
        print(f"SME response basis: {basis_cache_file(int_bins)}", flush=True)
        load_or_compute_basis(int_bins, jobs=args.jobs)
        print("SME response basis done")


if __name__ == "__main__":
    main()
