import argparse
import os

import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

from wremnants.postprocessing.datagroups.datagroups import Datagroups
from wremnants.production.muon_efficiencies_insitu import (
    insitu_n_coeff_pt,
    insitu_n_coeff_ut,
    insitu_pt_range,
    insitu_ut_range,
)
from wremnants.utilities import common
from wremnants.utilities.io_tools import base_io
from wums import logging, output_tools, plot_tools

hep.style.use(hep.style.ROOT)

# Produce the basis-orthogonalisation input for the in-situ muon efficiency
# method, from the same 4-category probe spectra that feed make_insitu_effMC.py.
#
# WHY. The coefficients are Chebyshev polynomials in (pT, uT), orthogonal under
# the analytic weight 1/sqrt(1-x^2) over the full x window. The probes are
# distributed nothing like that weight, so the per-block Fisher matrix is far
# from diagonal: postfit correlations reach |rho| ~ 0.99 and block condition
# numbers 1e4-1e5, which is what makes the fit need repeated refits to converge.
#
# WHAT. For each block -- (eta, charge) for idip/trigger, eta for the
# charge-inclusive iso -- we form the Gram matrix of the raw basis b under the
# per-bin weight that the fit actually sees,
#     G = sum_bins w_bin b b^T / sum_bins w_bin,
# and return A = L^-1 from the Cholesky factorisation G = L L^T. The
# transformed basis b' = A b then satisfies sum_bins w_bin b' b'^T = I, i.e. the
# per-block Fisher information (hence the postfit covariance) is ~diagonal.
#
# The weight is the squared first-order response of each probe population,
#     w = n_pass * 1 + n_fail * (eMC / (1 - eMC))^2,
# since at theta = 0 the log-derivative of a passing leg is b and that of a
# failing leg is -eMC*b/(1-eMC). eMC is recomputed here from the same spectra
# (the raw per-bin ratios of make_insitu_effMC.py); this only has to be
# approximately right -- see below.
#
# WHY AN APPROXIMATION IS FINE. Any invertible A spans the same space of
# scale-factor functions, so the fitted SF and its uncertainty band are
# invariant: A affects conditioning only, never the result. What it must be is
# FROZEN and SHARED -- the W and Z analyses correlate same-named coefficient
# nuisances, so they must use the same transform, which is why this is a
# standalone frozen input file rather than something each histmaker derives.
#
# Normalising G by the total weight keeps the coefficients on the same scale as
# the raw Chebyshev ones (G_00 = 1 since b_0 = 1, hence A_00 = 1 and the
# constant term is untouched). A is lower triangular, so b'_j mixes only
# b_0..b_j and the graded meaning of the coefficient index survives.
#
# ROW NORMALISATION. L^-1 alone is not enough. Its rows are large exactly along
# the weakly-constrained directions, so the orthonormalised b'_j -- which have
# unit RMS *under the density* -- can be huge out in the sparse tails, and the
# fail-category probes populate precisely those tails. The stored response is
# exp(delta * dlnW/dtheta), so that amplification pushes the log-normal
# linearisation far outside its validity (measured: max |delta*grad| 11 -> 61,
# i.e. exp() reaching 1e26). We therefore rescale each row so that
#     max_x |b'_j(x)| = 1 over the whole clamped x window,
# matching the raw Chebyshev bound (|T_k| <= 1), so the response magnitude is
# never worse than the current basis. Rescaling rows is a diagonal similarity
# on the Gram, so the Fisher stays exactly diagonal -- the correlations are
# still removed, only the (arbitrary, since unconstrained) nuisance units change.

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--inputFile",
    type=str,
    required=True,
    help="Histmaker output with the effMCprobe spectra (the --makeInsituEffMC "
    "pre-step output), used as the reference probe density",
)
parser.add_argument(
    "--process",
    type=str,
    default="Zmumu",
    help="Process whose probe spectra define the reference density. Use ONE "
    "reference for all analyses sharing the coefficient nuisances",
)
parser.add_argument(
    "--outpath", type=str, default=".", help="Output directory for the pkl"
)
parser.add_argument("--postfix", type=str, default="", help="Output filename postfix")
parser.add_argument(
    "--ridge",
    type=float,
    default=1e-10,
    help="Relative ridge added to the Gram diagonal before the Cholesky, to "
    "keep near-degenerate blocks factorisable",
)
parser.add_argument(
    "--plotdir", type=str, default=None, help="If set, write diagnostic plots here"
)
parser.add_argument("--debug", action="store_true", help="Print debug output")
args = parser.parse_args()

logger = logging.setup_logger(__file__, 4 if args.debug else 3)

CATEGORIES = ["nominal", "failIso", "failHLT", "failID"]
STEPS = ["idip", "trigger", "iso"]


def chebyshev(x, n):
    """T_0..T_{n-1} evaluated on x, shape (n, len(x))."""
    T = [np.ones_like(x), x]
    while len(T) < n:
        T.append(2 * x * T[-1] - T[-2])
    return np.array(T[:n])


def xtilde(v, lo, hi):
    return np.clip(2 * (np.asarray(v, dtype=float) - lo) / (hi - lo) - 1, -1, 1)


def read_probe_spectra(input_file, process):
    """Per-category probe counts (eta, pt, charge, genUT) plus the pt/uT axes."""
    datagroups = Datagroups(input_file)
    if process not in datagroups.groups or not datagroups.groups[process].members:
        raise ValueError(f"process '{process}' not available in {input_file}")
    for cat in CATEGORIES:
        datagroups.loadHistsForDatagroups(
            f"effMCprobe_{cat}", syst="", procsToRead=[process]
        )
    hists = {
        cat: datagroups.groups[process].hists[f"effMCprobe_{cat}"] for cat in CATEGORIES
    }
    names = tuple(a.name for a in hists["nominal"].axes)
    expected = ("eta", "pt", "charge", "genUT")
    if names != expected:
        raise ValueError(f"effMCprobe axes {names} != expected {expected}")
    n = {cat: hists[cat].values(flow=False) for cat in CATEGORIES}
    eta_ax, pt_ax, _, ut_ax = hists["nominal"].axes
    return n, eta_ax, pt_ax, ut_ax


def response_weights(n):
    """Per-bin weight w = n_pass + n_fail*(eMC/(1-eMC))^2 for each step, with
    the axes each step is differential in (idip has no uT, iso no charge).
    """

    def ratio(num, den):
        out = np.zeros_like(num, dtype=float)
        np.divide(num, den, out=out, where=den > 0)
        return out

    pass2HLT = n["nominal"] + n["failIso"]
    passID = pass2HLT + n["failHLT"]
    allProbes = passID + n["failID"]

    def fail_factor(eff):
        # (eMC/(1-eMC))^2, guarded against eMC -> 1 in empty/low-stat bins
        denom = np.clip(1.0 - eff, 1e-4, None)
        return (eff / denom) ** 2

    # idip: (eta, pt, charge), integrate uT away
    eff_idip = ratio(passID.sum(3), allProbes.sum(3))
    w_idip = passID.sum(3) + n["failID"].sum(3) * fail_factor(eff_idip)
    # trigger: (eta, pt, charge, uT)
    eff_trig = ratio(pass2HLT, passID)
    w_trig = pass2HLT + n["failHLT"] * fail_factor(eff_trig)
    # iso: (eta, pt, uT), charge inclusive
    eff_iso = ratio(n["nominal"].sum(2), pass2HLT.sum(2))
    w_iso = n["nominal"].sum(2) + n["failIso"].sum(2) * fail_factor(eff_iso)
    return {"idip": w_idip, "trigger": w_trig, "iso": w_iso}


def cholesky_transform(gram, ridge, tag):
    """A = L^-1 from gram = L L^T, with a relative ridge for stability.
    Falls back to the identity (no orthogonalisation) if the block is not
    factorisable at all, so a pathological bin never breaks the production.
    """
    n = gram.shape[0]
    scale = np.trace(gram) / n
    if not np.isfinite(scale) or scale <= 0:
        logger.warning(f"{tag}: empty/degenerate Gram matrix, using identity")
        return np.eye(n)
    g = gram + ridge * scale * np.eye(n)
    try:
        L = np.linalg.cholesky(g)
        return np.linalg.inv(L)
    except np.linalg.LinAlgError:
        logger.warning(f"{tag}: Cholesky failed, using identity")
        return np.eye(n)


def correlation_cond(gram):
    """Condition number of the Gram *correlation* matrix -- the basis-scale
    invariant measure of how degenerate the block is, and the one that governs
    how hard the block is to fit.
    """
    d = np.sqrt(np.diag(gram))
    if not np.all(d > 0):
        return np.nan
    return np.linalg.cond(gram / np.outer(d, d))


def row_normalise(A, basis_dense):
    """Scale each row of A so max_x |(A b)_j(x)| = 1 on the dense x grid, i.e.
    the same bound the raw Chebyshev basis has. Diagonal rescaling keeps the
    transformed Gram diagonal, so the decorrelation is untouched.
    """
    peak = np.abs(A @ basis_dense).max(axis=1)
    peak[~(peak > 0)] = 1.0
    return A / peak[:, None]


def build_transforms(w, pt_ax, ut_ax, n_eta, ridge):
    """Per-block A matrices, in the block order of the flat parameter layout
    (idip/trigger: charge-major then eta; iso: eta).
    """
    Tpt = chebyshev(xtilde(pt_ax.centers, *insitu_pt_range), insitu_n_coeff_pt)
    Tut = chebyshev(xtilde(ut_ax.centers, *insitu_ut_range), insitu_n_coeff_ut)
    # raw basis on the (pt) and (pt, uT) grids
    b_pt = Tpt  # (kPt, npt)
    b_2d = np.einsum("kp,mu->kmpu", Tpt, Tut).reshape(
        insitu_n_coeff_pt * insitu_n_coeff_ut, -1
    )  # (k2D, npt*nut)

    # dense x grids spanning the full clamped window, used only to bound the
    # transformed basis (the helper clamps to this window, so this is a true
    # bound on what any probe can produce)
    xd = np.linspace(-1.0, 1.0, 201)
    Td_pt = chebyshev(xd, insitu_n_coeff_pt)
    Td_ut = chebyshev(xd, insitu_n_coeff_ut)
    dense_pt = Td_pt
    dense_2d = np.einsum("kp,mu->kmpu", Td_pt, Td_ut).reshape(
        insitu_n_coeff_pt * insitu_n_coeff_ut, -1
    )

    n_2d = insitu_n_coeff_pt * insitu_n_coeff_ut
    out = {
        "idip": np.zeros((2 * n_eta, insitu_n_coeff_pt, insitu_n_coeff_pt)),
        "trigger": np.zeros((2 * n_eta, n_2d, n_2d)),
        "iso": np.zeros((n_eta, n_2d, n_2d)),
    }
    conds = {s: {"before": [], "after": []} for s in STEPS}

    def block(basis, dense, weights, tag):
        tot = weights.sum()
        if tot <= 0:
            logger.warning(f"{tag}: zero probe weight, using identity")
            return np.eye(basis.shape[0]), np.nan, np.nan
        gram = (basis * weights) @ basis.T / tot
        A = row_normalise(cholesky_transform(gram, ridge, tag), dense)
        after = A @ gram @ A.T
        # correlation, not covariance: the row scale is deliberate (see header)
        d = np.sqrt(np.diag(after))
        corr = after / np.outer(d, d)
        return A, correlation_cond(gram), np.linalg.cond(corr)

    for qbit in range(2):
        for b in range(n_eta):
            idx = qbit * n_eta + b
            A, c0, c1 = block(
                b_pt, dense_pt, w["idip"][b, :, qbit], f"idip eta{b} q{qbit}"
            )
            out["idip"][idx] = A
            conds["idip"]["before"].append(c0)
            conds["idip"]["after"].append(c1)

            A, c0, c1 = block(
                b_2d,
                dense_2d,
                w["trigger"][b, :, qbit, :].reshape(-1),
                f"trigger eta{b} q{qbit}",
            )
            out["trigger"][idx] = A
            conds["trigger"]["before"].append(c0)
            conds["trigger"]["after"].append(c1)

    for b in range(n_eta):
        A, c0, c1 = block(b_2d, dense_2d, w["iso"][b].reshape(-1), f"iso eta{b}")
        out["iso"][b] = A
        conds["iso"]["before"].append(c0)
        conds["iso"]["after"].append(c1)

    return out, conds


def plot_conditioning(plotdir, conds, n_eta):
    outdir = output_tools.make_plot_dir(plotdir)
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = {"idip": "tab:blue", "trigger": "tab:orange", "iso": "tab:red"}
    for step in STEPS:
        before = np.array(conds[step]["before"])
        after = np.array(conds[step]["after"])
        x = np.arange(len(before))
        ax.plot(x, before, color=colors[step], label=f"{step} (Chebyshev)")
        ax.plot(x, after, color=colors[step], ls="--", label=f"{step} (orthogonal)")
    ax.set_yscale("log")
    ax.set_xlabel("block index")
    ax.set_ylabel("Gram condition number")
    ax.legend(loc="upper right", ncol=2, fontsize="small")
    plot_tools.add_decor(ax, "CMS", "Simulation", data=False, lumi=None, loc=0)
    plot_tools.save_pdf_and_png(outdir, "insitu_basis_conditioning")
    output_tools.write_logfile(outdir, "insitu_basis_conditioning", args=args)
    logger.info(f"Wrote diagnostic plot to {outdir}")


def main():
    n, eta_ax, pt_ax, ut_ax = read_probe_spectra(args.inputFile, args.process)
    n_eta = eta_ax.size
    logger.info(
        f"Reference density from '{args.process}' in {args.inputFile}: "
        f"eta={n_eta}, pt={pt_ax.size}, uT={ut_ax.size}"
    )

    w = response_weights(n)
    transforms, conds = build_transforms(w, pt_ax, ut_ax, n_eta, args.ridge)

    for step in STEPS:
        before = np.array(conds[step]["before"])
        after = np.array(conds[step]["after"])
        logger.info(
            f"{step}: Gram correlation condition number median {np.nanmedian(before):.1f} "
            f"-> {np.nanmedian(after):.3f} (max {np.nanmax(before):.1f} "
            f"-> {np.nanmax(after):.3f}) over {len(before)} blocks"
        )

    if args.plotdir:
        plot_conditioning(args.plotdir, conds, n_eta)

    outname = "insitu_basis" + (f"_{args.postfix}" if args.postfix else "")
    outfile = os.path.join(args.outpath, f"{outname}.pkl.lz4")
    output_dict = {
        "transform": transforms,
        "reference": args.process,
        "n_eta": n_eta,
        "n_coeff_pt": insitu_n_coeff_pt,
        "n_coeff_ut": insitu_n_coeff_ut,
        "pt_range": tuple(insitu_pt_range),
        "ut_range": tuple(insitu_ut_range),
    }
    output_tools.write_lz4_pkl_output(
        outfile,
        outname,
        output_dict,
        common.base_dir,
        args,
        file_meta_data=base_io.get_metadata(args.inputFile),
    )
    logger.info(f"Wrote in-situ basis transform to {outfile}")


if __name__ == "__main__":
    main()
