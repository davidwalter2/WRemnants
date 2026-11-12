#!/usr/bin/env python3
"""Expected 1D bounds on the sidereal d=4 SME coefficients from the ATLAS 13 TeV Z-pole
cross section, as in Table 4 of Lunghi & Sherrill, arXiv:2412.14305, with rabbit.

Observable: dsigma/dm_ll evaluated exactly at m_ll = m_Z (not integrated over an m_ll bin),
from the double differential SME model (sme_double_differential.py, LIV_mll_pheno),
integrated over the dilepton rapidity:
    --yllMax all   all rapidities (kinematic limit |y| < ln(sqrt(s) / m_Z))
    --yllMax 2.5   |y_ll| < 2.5
The value is the same for each of --nTimeBins sidereal time bins (default 4, as the paper);
the LIV contribution delta(dsigma/dm_ll)(t), averaged over each time bin, is the sidereal
modulation. The pseudo-data are the SM prediction (asimov), so the fitted coefficients are 0
and their uncertainties are the expected bounds: the coefficient for which the additional
LIV contribution to the cross section is as large as the ATLAS uncertainty.

The model is told to evaluate at a point through the 'liv_mll_point' meta data (the channel
has only a time axis); the basis (basis.bin_xsec) then computes 2 m dsigma/dm^2 at m = m_Z
instead of the integral over a bin. Equivalent by hand: --paramModel LIV_mll_pheno mllPoint=91.1879.

Uncertainties: the ATLAS measurement behind the paper, arXiv:1912.02844 (Eq. 16 of 2412.14305),
    sigma_exp  = 736.2 +- 0.2 (stat) +- 6.4 (syst) +- 15.5 (lumi) pb
and the NNLO prediction of arXiv:1612.03636 (Eq. 17), used for the theory uncertainty,
    sigma_NNLO = 703 +19/-24 +6/-8 +4/-6 +5/-5 pb.
They enter as relative uncertainties on dsigma/dm_ll in each time bin:
  - stat: uncorrelated between the time bins, scaled by sqrt(nTimeBins) (each bin has
    1/nTimeBins of the events); --uncorrScaling none (or --noStatScaling) keeps the unscaled
    0.027%, --uncorrScaling all scales the uncorrelated systematics by sqrt(nTimeBins) too. This is the
    diagonal of the chi2. The fit must use --covarianceFit, which takes the covariance from
    these data variances; --chisqFit in rabbit_fit.py uses the yields as variances (set_nobs
    is called without variances for asimov and data) and would ignore the ATLAS uncertainty.
  - th, lumi, syst: by default (--correlated th lumi syst, the last column of Table 4) 100%
    correlated between the time bins, one normalisation nuisance per source acting on all
    time bins. Sources not listed in --correlated are added in quadrature to the per-bin
    (uncorrelated) uncertainty, without the sqrt(n) scaling, giving the other columns:
        --correlated th           -> column "delta_th"
        --correlated th lumi      -> column "delta_th, delta_lumi"
        --correlated th lumi syst -> column "delta_th, delta_lumi, delta_syst" (default)

Comparison with Table 4: the sidereal shifts sum to zero over the time bins, so the correlated
sources drop out and every bound is proportional to the per-bin uncorrelated uncertainty.
The paper's column ratios (delta_th,lumi / all correlated = 32 = sqrt(stat^2 + syst^2) / stat
with the unscaled stat) need stat and the uncorrelated syst scaled alike: --uncorrScaling none
or all (the two differ by an overall factor 2 only). The default (stat only, as the paper's
text says) gives a ratio 16. The paper's bounds are, as in arXiv:1610.08755 / 1805.11684, the
median over pseudo-experiments of the 95% CL upper bound (|c_hat| + 1.96 sigma = 2.63 sigma
median, 2.76 sigma mean). With --fit, --nToys pseudo-experiments are fitted individually with
rabbit (--toys, normal data randomisation from the covariance, frequentist randomisation of the
correlated nuisances) and the asimov 95% CI, the mean and the median bound are printed.

With --fit the script also prints three cross-checks computed directly from the model
templates, without rabbit's minimiser:
  - chi2: sigma_c = (s^T V^-1 s)^-1/2, s the relative LIV shift per time bin for a unit
    coefficient, V the relative covariance (uncorrelated + correlated sources); this should
    equal the rabbit result.
  - sat: the saturation bound delta_bin / max_t |s_t|, i.e. the LIV contribution in the most
    sensitive time bin equal to the per-bin uncorrelated uncertainty.
  - ana. mean / ana. med: the mean / median 95% bound of analytic linear pseudo-experiments.

  usage (inside singularity, after source setup.sh):
    python scripts/plotting/liv_zpole_atlas_xsec_fit.py -o OUTDIR --fit all
    python scripts/plotting/liv_zpole_atlas_xsec_fit.py -o OUTDIR --yllMax 2.5 --fit all
    python scripts/plotting/liv_zpole_atlas_xsec_fit.py -o OUTDIR --correlated th lumi --fit all
    # a single fit by hand
    rabbit_fit.py OUTDIR/liv_zpole_atlas_mZ_yllall_4timeBins_corr_th_lumi_syst.hdf5 \\
        --paramModel LIV_mll_pheno cxyu --covarianceFit --noBinByBinStat -o OUTDIR
"""

import argparse
import os
import subprocess

import hist
import numpy as np

from rabbit import tensorwriter
from rabbit.param_models.liv_rapidity.basis import YLL_FULL, load_or_compute_basis
from rabbit.param_models.liv_rapidity.functions import m_Z
from rabbit.param_models.sme_double_differential import SU2_ALIASES, SU2_BASIS

# arXiv:1912.02844 (Eq. 16 of arXiv:2412.14305), fiducial Z -> ll cross section [pb]
ATLAS_XSEC = 736.2
ATLAS_STAT = 0.2
ATLAS_SYST = 6.4
ATLAS_LUMI = 15.5
# arXiv:1612.03636 (Eq. 17 of arXiv:2412.14305), NNLO prediction [pb], (up, down) per component
NNLO_XSEC = 703.0
NNLO_UNC = [(19.0, 24.0), (6.0, 8.0), (4.0, 6.0), (5.0, 5.0)]

SOURCES = ["th", "lumi", "syst"]
GEV2_TO_PB = 0.3893794e9  # (hbar c)^2 in pb GeV^2, the basis is in natural units
DIRECTIONS = ["xy", "xz", "yz", "xx"]  # order of Table 4
SCALE = 1e-6  # LIV_rapidity parameter scale


def all_coefficients():
    """Table 4 coefficients available in the SU(2) basis of the model"""
    return [
        f"{mode}{d}{q}"
        for mode, q in SU2_BASIS
        if (mode, q) not in SU2_ALIASES
        for d in DIRECTIONS
    ]


def label(c):
    d = "XX-YY" if c[1:3] == "xx" else c[1:3].upper()
    return f"|{c[0]}_{c[3]}^{d}|"


def uncertainties(args):
    """relative uncertainties: per-bin uncorrelated, {correlated source: value}, theory components"""
    n = args.nTimeBins
    th_comps = [(u / NNLO_XSEC, d / NNLO_XSEC) for u, d in NNLO_UNC]
    rel = {
        "syst": ATLAS_SYST / ATLAS_XSEC,
        "lumi": ATLAS_LUMI / ATLAS_XSEC,
    }
    if not args.noTheory:
        rel["th"] = np.sqrt(sum((0.5 * (u + d)) ** 2 for u, d in th_comps))
    scale = np.sqrt(n)
    stat = ATLAS_STAT / ATLAS_XSEC * (1.0 if args.uncorrScaling == "none" else scale)
    correlated = {s: rel[s] for s in args.correlated if s in rel}
    uncorr_syst2 = sum(v**2 for s, v in rel.items() if s not in correlated)
    if args.uncorrScaling == "all":
        uncorr_syst2 *= scale**2
    rel_bin = np.sqrt(stat**2 + uncorr_syst2)
    return stat, rel, correlated, rel_bin, th_comps


def write_input(args):
    n = args.nTimeBins
    ymax = YLL_FULL if args.yllMax.lower() == "all" else float(args.yllMax)
    m = float(args.mllPoint)

    # SM dsigma/dm_ll at m_ll = m (point bin), from the same basis the model will use
    sm_total, _ = load_or_compute_basis([(m, m, -ymax, ymax)], jobs=args.jobs)
    dsig = float(sm_total[0]) * GEV2_TO_PB  # [pb/GeV]

    stat, rel, correlated, rel_bin, th_comps = uncertainties(args)
    print(
        f"dsigma/dm_ll at m_ll = {m:g} GeV, |y_ll| < {ymax:.4g}, {n} sidereal time bins"
    )
    print(f"  SM (LO, model) = {dsig:.6g} pb/GeV in every time bin")
    print(
        f"  stat : {stat:.3e}  (uncorrelated{'' if args.uncorrScaling == 'none' else f', x sqrt({n})'})"
    )
    for s, v in rel.items():
        if s in correlated:
            print(f"  {s:5s}: {v:.3e}  (correlated)")
        else:
            sc = (
                f", x sqrt({n}) in the per-bin sum"
                if args.uncorrScaling == "all"
                else ""
            )
            print(f"  {s:5s}: {v:.3e}  (uncorrelated{sc})")
    print(f"  per-bin uncorrelated relative uncertainty: {rel_bin:.3e}")

    axis_time = hist.axis.Regular(
        n, 0, 24, name="time", overflow=False, underflow=False
    )
    h = hist.Hist(axis_time, storage=hist.storage.Weight())
    h.values()[...] = dsig
    h.variances()[...] = 0.0  # no MC stat uncertainty on the prediction

    writer = tensorwriter.TensorWriter(systematic_type="log_normal")
    ch = "ch_z"
    writer.add_channel(h.axes, ch)
    data_var = np.full(n, (rel_bin * dsig) ** 2)
    writer.add_data(h.values(), ch, variances=data_var)
    # explicit (diagonal) covariance for --covarianceFit: the uncorrelated per-bin part only,
    # the correlated sources are nuisance parameters below
    writer.add_data_covariance(np.diag(data_var))
    writer.add_process(h, "Zll", ch, signal=True)

    # 100% correlated between the time bins: one normalisation nuisance per source
    for s in correlated:
        if s == "th":
            for i, (u, d) in enumerate(th_comps):
                writer.add_norm_systematic(
                    f"atlas_th{i}", "Zll", ch, [(1 + u, 1 - d)], groups=["th"]
                )
        else:
            writer.add_norm_systematic(f"atlas_{s}", "Zll", ch, 1 + rel[s], groups=[s])

    yll_tag = "all" if args.yllMax.lower() == "all" else f"{ymax:g}".replace(".", "p")
    corr_tag = "_".join(correlated) if correlated else "none"
    outname = (
        args.outname or f"liv_zpole_atlas_mZ_yll{yll_tag}_{n}timeBins_corr_{corr_tag}"
    )
    if args.noTheory:
        outname += "_noTh"
    if args.uncorrScaling != "stat":
        outname += f"_uncorrScaling{args.uncorrScaling.capitalize()}"
    if args.postfix:
        outname += f"_{args.postfix}"
    writer.write(
        outfolder=args.outdir,
        outfilename=outname,
        meta_data_dict={
            "liv_experiment": args.experiment,
            "liv_yll_max": ymax,
            "liv_mll_point": m,
        },
    )
    path = os.path.join(args.outdir, outname + ".hdf5")
    print(f"wrote {path}")
    return path, outname


def direct_bounds(args, path, coeffs):
    """chi2 and saturation bounds from the model templates (no minimisation)"""
    from rabbit import inputdata
    from rabbit.param_models import helpers as ph

    indata = inputdata.FitInputData(path)
    model = ph.load_model("LIV_mll_pheno", indata, *coeffs)
    names = [p.decode() for p in model.params]
    # relative shift per time bin for a unit coefficient [npoi, ntime]
    # (model.sme is the change in cross section per unit coefficient, the parameter is
    # coefficient / SCALE)
    shift = model.sme.numpy() / model.sm_sigma.numpy()[None]
    _, _, correlated, rel_bin, _ = uncertainties(args)
    n = shift.shape[1]
    V = np.diag(np.full(n, rel_bin**2)) + sum(
        v**2 for v in correlated.values()
    ) * np.ones((n, n))
    Vinv = np.linalg.inv(V)
    rng = np.random.default_rng(args.seed)
    # pseudo-data: relative fluctuations of the time bins with covariance V
    toys = rng.multivariate_normal(np.zeros(n), V, size=args.nToys)
    out = {}
    for c in coeffs:
        s = shift[names.index(c)]
        sigma = 1.0 / np.sqrt(s @ Vinv @ s)
        # per pseudo-experiment: chi2 fit c_hat = sigma^2 s^T V^-1 d, 95% CL interval
        # c_hat +- 1.96 sigma, upper bound on |c| = |c_hat| + 1.96 sigma; take the median
        chat = sigma**2 * (toys @ Vinv @ s)
        bound = np.abs(chat) + 1.96 * sigma
        out[c] = (
            sigma,
            rel_bin / np.max(np.abs(s)),
            s,
            np.median(bound),
            np.mean(bound),
        )
    return out


def read_toys(filename, poi, ntoys):
    """fitted value and Hessian uncertainty of poi in each pseudo-experiment [coefficient units]"""
    import h5py

    from rabbit import io_tools

    vals, errs = [], []
    with h5py.File(filename, mode="r") as f:
        for itoy in range(1, ntoys + 1):
            parms = io_tools.get_fitresult(f, f"toy{itoy}")["parms"].get()
            i = list(parms.axes[0]).index(poi)
            vals.append(parms.values()[i])
            errs.append(np.sqrt(parms.variances()[i]))
    return np.array(vals) * SCALE, np.array(errs) * SCALE


def run_fits(args, path, outname):
    from rabbit import io_tools

    coeffs = all_coefficients() if args.fit == ["all"] else args.fit
    fitdir = os.path.join(args.outdir, f"fits_{outname}")
    os.makedirs(fitdir, exist_ok=True)
    ntoys = args.nToys if not args.noToyFits else 0
    results, toys = {}, {}
    for c in coeffs:
        fname = f"fitresults_{c}.hdf5"
        cmd = [
            "rabbit_fit.py",
            path,
            "--paramModel",
            "LIV_mll_pheno",
            c,
            "--covarianceFit",
            "--noBinByBinStat",
            "-o",
            fitdir,
            "--outname",
            fname,
        ]
        if ntoys > 0:
            # asimov fit plus ntoys pseudo-experiments, each fitted on its own: the data are drawn
            # from the uncorrelated covariance (normal, as required for --covarianceFit) and the
            # constraint minima of the correlated nuisances are randomised (frequentist)
            cmd += [
                "--toys",
                "-1",
                str(ntoys),
                "--toysDataRandomize",
                "normal",
                "--toysSystRandomize",
                "frequentist",
                "--seed",
                str(args.seed),
            ]
        cmd += args.fitArgs
        print(" ".join(cmd), flush=True)
        out = subprocess.run(cmd, capture_output=True, text=True)
        if out.returncode != 0:
            print(out.stdout[-3000:], out.stderr[-3000:])
            raise SystemExit(f"fit for {c} failed")
        fpath = os.path.join(fitdir, fname)
        parms = io_tools.get_fitresult(fpath, "asimov" if ntoys > 0 else None)[
            "parms"
        ].get()
        i = list(parms.axes[0]).index(c)
        results[c] = (parms.values()[i] * SCALE, np.sqrt(parms.variances()[i]) * SCALE)
        if ntoys > 0:
            toys[c] = read_toys(fpath, c, ntoys)

    direct = direct_bounds(args, path, coeffs)

    print(
        f"\nexpected bounds, dsigma/dm_ll at m_ll = {args.mllPoint:g} GeV, |y_ll| < {args.yllMax}, "
        f"correlated: {args.correlated}, uncorrScaling: {args.uncorrScaling}, {args.nTimeBins} time bins"
    )
    if ntoys > 0:
        print(
            f"  95% CI: asimov fit, c_hat = 0 +- 1.96 sigma.  For each of {ntoys} pseudo-experiments, "
            "fitted individually with rabbit, the 95% CL bound on |c| is |c_hat| + 1.96 sigma;"
        )
        print(
            "  +mean / +median: mean / median of that bound over the pseudo-experiments "
            "(Gaussian expectation (0.798 + 1.96) sigma / (0.674 + 1.96) sigma, in brackets);"
        )
        print(
            "  toy rms/sig: spread of c_hat over the pseudo-experiments / asimov sigma (should be 1)"
        )
        print(
            f"{'coefficient':13s} {'95% CI':>9s} {'+mean':>9s} {'+median':>9s}"
            f" {'(+mean)':>9s} {'(+median)':>9s} {'toy rms/sig':>11s}"
        )
        for c, (v, e) in results.items():
            chat, chat_err = toys[c]
            bound = np.abs(chat) + 1.96 * chat_err
            print(
                f"{label(c):13s} {1.96 * e:9.2e} {np.mean(bound):9.2e} {np.median(bound):9.2e}"
                f" {(0.7979 + 1.96) * e:9.2e} {(0.6745 + 1.96) * e:9.2e} {np.std(chat) / e:11.3f}"
            )
        print()

    print(
        "  cross-checks: rabbit 1 sigma and 95% CL (1.96 sigma) of the asimov fit; chi2: analytic sigma; "
        "sat: LIV shift in the most sensitive time bin = per-bin uncertainty;"
    )
    print(
        f"  ana. mean / ana. med: the same mean / median bounds from {args.nToys} analytic linear "
        "pseudo-experiments (no rabbit)"
    )
    print(
        f"{'coefficient':13s} {'rabbit 1sig':>11s} {'rabbit 95%':>11s} {'chi2 1sig':>11s} {'sat':>11s}"
        f" {'ana. mean':>11s} {'ana. med':>11s}   rel. shift per time bin for c = 1e-4"
    )
    for c, (v, e) in results.items():
        chi2, sat, s, med, mean = direct[c]
        shifts = " ".join(f"{x * 1e-4:+.2e}" for x in s)
        print(
            f"{label(c):13s} {e:11.2e} {1.96 * e:11.2e} {chi2:11.2e} {sat:11.2e} {mean:11.2e}"
            f" {med:11.2e}   {shifts}"
        )
    save = dict(
        coefficients=np.array(list(results)),
        values=np.array([v for v, _ in results.values()]),
        sigma=np.array([e for _, e in results.values()]),
        sigma_chi2=np.array([direct[c][0] for c in results]),
        saturation=np.array([direct[c][1] for c in results]),
        shift_per_unit=np.array([direct[c][2] for c in results]),
        median_analytic_bound95=np.array([direct[c][3] for c in results]),
        mean_analytic_bound95=np.array([direct[c][4] for c in results]),
    )
    if ntoys > 0:
        save["toy_values"] = np.array([toys[c][0] for c in results])
        save["toy_sigma"] = np.array([toys[c][1] for c in results])
    np.savez(os.path.join(fitdir, "bounds.npz"), **save)
    print(f"saved {os.path.join(fitdir, 'bounds.npz')}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("-o", "--outdir", default=os.getcwd())
    ap.add_argument("--outname", default=None)
    ap.add_argument("-p", "--postfix", default=None)
    ap.add_argument("--nTimeBins", type=int, default=4)
    ap.add_argument(
        "--mllPoint",
        type=float,
        default=m_Z,
        help="m_ll [GeV] of dsigma/dm_ll (default: m_Z)",
    )
    ap.add_argument(
        "--yllMax",
        default="all",
        help="integrate |y_ll| < this: 'all' (full rapidity range, default) or a number, e.g. 2.5",
    )
    ap.add_argument(
        "--correlated",
        nargs="*",
        default=SOURCES,
        choices=SOURCES,
        help="sources taken 100%% correlated between time bins (default: all, last column of Table 4)",
    )
    ap.add_argument(
        "--noTheory", action="store_true", help="drop the NNLO theory uncertainty"
    )
    ap.add_argument(
        "--uncorrScaling",
        default="stat",
        choices=["stat", "all", "none"],
        help="per-bin uncorrelated uncertainties scaled by sqrt(nTimeBins): 'stat' only (default), "
        "'all' (stat and the uncorrelated systematics), 'none'. Table 4 of the paper is reproduced "
        "by 'all' or 'none' (only the stat/syst ratio matters; see the docstring)",
    )
    ap.add_argument(
        "--noStatScaling",
        action="store_const",
        dest="uncorrScaling",
        const="none",
        help="same as --uncorrScaling none",
    )
    ap.add_argument(
        "--nToys",
        type=int,
        default=1000,
        help="pseudo-experiments, each fitted with rabbit, for the mean / median 95%% CL bounds",
    )
    ap.add_argument(
        "--noToyFits",
        action="store_true",
        help="skip the rabbit pseudo-experiment fits (only the asimov fit and the analytic toys)",
    )
    ap.add_argument(
        "--seed", type=int, default=12345, help="random seed of the pseudo-experiments"
    )
    ap.add_argument(
        "--experiment", default="ATLAS", type=str.upper, choices=["CMS", "ATLAS"]
    )
    ap.add_argument(
        "--jobs", type=int, default=8, help="workers for the SME basis if not cached"
    )
    ap.add_argument(
        "--fit",
        nargs="*",
        default=None,
        help="run one 1D rabbit fit per coefficient ('all' for the Table 4 set) and print the bounds",
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
            f"fit with: rabbit_fit.py {path} --paramModel LIV_mll_pheno cxyu --covarianceFit --noBinByBinStat"
        )


if __name__ == "__main__":
    main()
