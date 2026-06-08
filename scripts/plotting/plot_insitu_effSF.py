import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import chebyshev as C

from rabbit.io_tools import get_fitresult
from wremnants.production.muon_efficiencies_insitu import (
    insitu_delta,
    insitu_effMC_max,
    insitu_n_coeff_pt,
    insitu_n_coeff_ut,
    insitu_parameter_labels,
    insitu_ut_range,
    load_insitu_central,
    make_muon_insitu_effMC_helper,
)
from wums import logging, output_tools, plot_tools

# Postfit diagnostic for the in-situ tag-and-probe muon efficiency fit.
#
# The fit floats the ID/HLT/Iso data/MC scale factors as Chebyshev polynomials
#   SF_X(pt[,uT]) = exp( sum_{k,m} theta_{X,eta,[q],k,m} T_k(x_pt) T_m(x_uT) ),
# with theta = theta_central + delta * n_hat (n_hat = fitted nuisance). The
# *marginal* per-coefficient uncertainties are large and strongly correlated
# within each (eta, charge) block, so they are not the meaningful quantity --
# the SF (a fixed linear combination of the coefficients) is what the data
# constrains. This script propagates the FULL parameter covariance into the SF
# band sigma_SF(pt[,uT]) and overlays the external tag-and-probe SF
# (= EffData2D / EffMC2D), in the same per-(eta, charge) pt-cell layout as the
# effMC diagnostic plots (make_insitu_effMC.py).

logger = logging.child_logger(__name__)

# Chebyshev coefficient block per step: prefix + whether it is charge-dependent
# and has a uT dependence. Mirrors insitu_parameter_labels().
STEP_INFO = {
    "idip": dict(prefix="effInsituID", charge=True, has_ut=False, ylabel="ID"),
    "trigger": dict(prefix="effInsituHLT", charge=True, has_ut=True, ylabel="trigger"),
    "iso": dict(prefix="effInsituIso", charge=False, has_ut=True, ylabel="iso"),
}
# step -> external T&P file token; idip/trigger charge-split, iso "both".
TNP_FILE = {"idip": "idip", "trigger": "trigger", "iso": "iso"}
TNP_STYLE = dict(fmt="s", color="black", markersize=4, capsize=2, zorder=5)
CHARGE_TAGS = {-1: ("minus", r"$\mu^{-}$"), +1: ("plus", r"$\mu^{+}$")}


def cheb(x, n):
    """Chebyshev basis [T_0(x), ..., T_{n-1}(x)] at scalar/array x."""
    return np.array([C.chebval(x, [0] * i + [1]) for i in range(n)])


def norm(v, lo, hi):
    return 2.0 * (np.clip(v, lo, hi) - lo) / (hi - lo) - 1.0


_tnp_cache = {}


def load_tnp_sf(tnp_dir, step, charge_tag):
    """External T&P scale factor SF = EffData2D / EffMC2D for one (step, charge).

    Returns (pt_centers[npt], pt_halfwidths[npt], sf[n_eta, npt], err[n_eta, npt])
    with eta on the first axis (the 48-bin grid, identical to ours) and pt on
    the second, or None if the file is missing. Cached per (step, charge).
    """
    import ROOT

    key = (tnp_dir, step, charge_tag)
    if key in _tnp_cache:
        return _tnp_cache[key]
    path = os.path.join(
        tnp_dir, f"allEfficiencies_2D_{TNP_FILE[step]}_{charge_tag}.root"
    )
    if not os.path.isfile(path):
        logger.warning(f"T&P file not found, skipping overlay: {path}")
        _tnp_cache[key] = None
        return None
    f = ROOT.TFile.Open(path)
    eff = {}
    for which, hname in (("mc", "EffMC2D"), ("data", "EffData2D")):
        h2 = f.Get(hname)
        n_eta, n_pt = h2.GetNbinsX(), h2.GetNbinsY()
        pt_c = np.array([h2.GetYaxis().GetBinCenter(j + 1) for j in range(n_pt)])
        pt_hw = np.array([0.5 * h2.GetYaxis().GetBinWidth(j + 1) for j in range(n_pt)])
        e = np.array(
            [
                [h2.GetBinContent(i + 1, j + 1) for j in range(n_pt)]
                for i in range(n_eta)
            ]
        )
        u = np.array(
            [[h2.GetBinError(i + 1, j + 1) for j in range(n_pt)] for i in range(n_eta)]
        )
        eff[which] = (pt_c, pt_hw, e, u)
    f.Close()
    pt_c, pt_hw, ed, errd = eff["data"]
    _, _, em, errm = eff["mc"]
    with np.errstate(divide="ignore", invalid="ignore"):
        sf = ed / em
        rel = np.sqrt((errd / ed) ** 2 + (errm / em) ** 2)
        err = np.abs(sf) * rel
    sf = np.where(np.isfinite(sf), sf, np.nan)
    err = np.where(np.isfinite(err), err, 0.0)
    out = (pt_c, pt_hw, sf, err)
    _tnp_cache[key] = out
    return out


def build_blocks(fitresult, prev_file, pt_range):
    """Per (step, eta, charge) Chebyshev coefficients theta and covariance.

    Returns (n_eta, blocks) where blocks[step] is a dict keyed by (eta, charge
    or None) -> (theta[ncoeff], cov[ncoeff, ncoeff]). theta already includes the
    previous central value; cov = delta^2 * cov(n_hat) (the fixed central offset
    does not contribute to the uncertainty).
    """
    parms = fitresult["parms"].get()
    labels = np.array(parms.axes["parms"]).astype(str)
    nhat = parms.values()
    cov = fitresult["cov"].get().values()
    idx = {l: i for i, l in enumerate(labels)}

    n_eta = 1 + max(
        int(l.split("_eta")[1].split("_")[0])
        for l in labels
        if l.startswith("effInsituIso_eta")
    )
    prev = dict(
        zip(
            insitu_parameter_labels(n_eta, insitu_n_coeff_pt, insitu_n_coeff_ut),
            load_insitu_central(prev_file, n_eta, insitu_n_coeff_pt, insitu_n_coeff_ut),
        )
    )

    blocks = {}
    for step, info in STEP_INFO.items():
        pref, has_ut = info["prefix"], info["has_ut"]
        charges = ("minus", "plus") if info["charge"] else (None,)
        nk, nm = insitu_n_coeff_pt, (insitu_n_coeff_ut if has_ut else 1)
        blocks[step] = {}
        for b in range(n_eta):
            for q in charges:
                names = []
                for k in range(nk):
                    for m in range(nm):
                        n = f"{pref}_eta{b}"
                        if q is not None:
                            n += f"_q{q}"
                        n += f"_cPt{k}"
                        if has_ut:
                            n += f"_cUt{m}"
                        names.append(n)
                ids = np.array([idx[n] for n in names])
                theta = np.array([prev[n] for n in names]) + insitu_delta * nhat[ids]
                cblock = insitu_delta**2 * cov[np.ix_(ids, ids)]
                blocks[step][(b, q)] = (theta, cblock)
    return n_eta, blocks


def effmc_along_pt(effmc, step, i_eta, qkey, ut, pt_grid):
    """MC efficiency e at (i_eta, pt_grid, [charge], [ut]) for one cell. Axes are
    indexed positionally: idip (eta,pt,charge); trigger (eta,pt,charge,ut);
    iso (eta,pt,ut)."""
    h = effmc[step]
    v = h.values()
    axes = list(h.axes)
    pt_ax = axes[1]
    pt_i = np.clip([pt_ax.index(float(pt)) for pt in pt_grid], 0, pt_ax.size - 1)
    if step == "idip":
        return v[i_eta, pt_i, axes[2].index(float(qkey))]
    if step == "trigger":
        return v[i_eta, pt_i, axes[2].index(float(qkey)), axes[3].index(float(ut))]
    return v[i_eta, pt_i, axes[2].index(float(ut))]  # iso


def sf_band(theta, cblock, pt_grid, ut, pt_range, ut_range, has_ut, e):
    """data/MC SF(pt) = eps/e and sigma_SF(pt) at fixed uT (logit param, full cov).

    eps = e*e^s / (e*e^s + 1-e) is the data efficiency; SF = eps/e = e^s/D.
    The covariance Jacobian carries d ln SF/d theta_c = (1-eps) * B_c, so
    sigma_lnSF = (1-eps) * sqrt(B^T Cov B). ``e`` is the effMC over pt_grid.
    """
    xpt = norm(pt_grid, *pt_range)
    Tp = cheb(xpt, insitu_n_coeff_pt)  # (nk, npt)
    if has_ut:
        Tu = cheb(norm(ut, *ut_range), insitu_n_coeff_ut)  # (nm,)
        B = (Tp[:, None, :] * Tu[None, :, None]).reshape(-1, pt_grid.size)  # (nc, npt)
    else:
        B = Tp
    expS = np.exp(theta @ B)
    e = np.clip(np.asarray(e, float), 1e-6, insitu_effMC_max)
    D = e * expS + (1.0 - e)
    eps = e * expS / D
    sf = expS / D  # = eps / e
    var = (1.0 - eps) ** 2 * np.einsum("ip,ij,jp->p", B, cblock, B)
    return sf, sf * np.sqrt(np.maximum(var, 0.0))


def make_plots(fitresult, args):
    n_eta, blocks = build_blocks(fitresult, args.prevSFFile, args.ptRange)
    effmc = make_muon_insitu_effMC_helper(args.effMCFile)
    outdir = output_tools.make_plot_dir(*args.plotdir.rsplit("/", 1), eoscp=args.eoscp)
    eta_edges = np.linspace(-2.4, 2.4, n_eta + 1)
    pt_grid = np.linspace(args.ptRange[0], args.ptRange[1], 120)
    ut_slices = args.utSlices
    ut_cmap = plt.cm.viridis(np.linspace(0.15, 0.85, len(ut_slices)))

    for step, info in STEP_INFO.items():
        has_ut = info["has_ut"]
        step_dir = os.path.join(outdir, step, "per_cell")
        os.makedirs(step_dir, exist_ok=True)
        charges = (-1, +1) if info["charge"] else (None,)
        for i_eta in range(n_eta):
            eta_lo, eta_hi = eta_edges[i_eta], eta_edges[i_eta + 1]
            eta_text = f"${eta_lo:+.2f}\\,\\leq\\,\\eta\\,<\\,{eta_hi:+.2f}$"
            for qkey in charges:
                q = None if qkey is None else CHARGE_TAGS[qkey][0]
                theta, cblock = blocks[step][(i_eta, q)]
                fig, ax = plt.subplots(figsize=(8, 6))

                if has_ut:
                    for ut, col in zip(ut_slices, ut_cmap):
                        e = effmc_along_pt(effmc, step, i_eta, qkey, ut, pt_grid)
                        sf, sig = sf_band(
                            theta,
                            cblock,
                            pt_grid,
                            ut,
                            args.ptRange,
                            insitu_ut_range,
                            True,
                            e,
                        )
                        ax.plot(pt_grid, sf, color=col, label=f"$u_T={ut:g}$ GeV")
                        ax.fill_between(
                            pt_grid, sf - sig, sf + sig, color=col, alpha=0.25, lw=0
                        )
                else:
                    e = effmc_along_pt(effmc, step, i_eta, qkey, 0.0, pt_grid)
                    sf, sig = sf_band(
                        theta,
                        cblock,
                        pt_grid,
                        0.0,
                        args.ptRange,
                        insitu_ut_range,
                        False,
                        e,
                    )
                    ax.plot(pt_grid, sf, color="#5790FC", label="in-situ SF")
                    ax.fill_between(
                        pt_grid,
                        sf - sig,
                        sf + sig,
                        color="#5790FC",
                        alpha=0.3,
                        lw=0,
                        label=r"$\pm\,\sigma_{SF}$ (full cov)",
                    )

                # external T&P SF overlay (uT-integrated 2D), reference points
                if args.tagAndProbeDir:
                    tag = "both" if q is None else q
                    tnp = load_tnp_sf(args.tagAndProbeDir, step, tag)
                    if tnp:
                        pt_c, pt_hw, sf_t, err_t = tnp
                        lbl = "T&P SF" + (" ($u_T$ incl.)" if has_ut else "")
                        ax.errorbar(
                            pt_c,
                            sf_t[i_eta],
                            yerr=err_t[i_eta],
                            xerr=pt_hw,
                            label=lbl,
                            **TNP_STYLE,
                        )

                annot = eta_text if q is None else f"{eta_text}\n{CHARGE_TAGS[qkey][1]}"
                suffix = "" if q is None else f"_q{q}"
                ax.text(
                    0.04,
                    0.06,
                    annot,
                    transform=ax.transAxes,
                    ha="left",
                    va="bottom",
                    fontsize=18,
                )
                ax.axhline(1.0, ls="--", color="gray", lw=0.8, zorder=0)
                ax.set_xlim(*args.ptRange)
                ax.set_ylim(*args.ylim)
                ax.set_xlabel(r"$p_T^\mu$ [GeV]")
                ax.set_ylabel(f"data/MC SF ({info['ylabel']})")
                ax.legend(loc="upper right", frameon=False, fontsize=12, ncol=2)
                plot_tools.add_decor(
                    ax, "CMS", "Preliminary", data=True, lumi=16.8, loc=0
                )
                name = f"effSF_{step}_pt_eta{i_eta:02d}{suffix}"
                plot_tools.save_pdf_and_png(step_dir, name)
                plt.close(fig)
                output_tools.write_index_and_log(step_dir, name, args=args)
        logger.info(f"Wrote {step} SF band plots to {step_dir}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "-i",
        "--inputFile",
        required=True,
        help="rabbit fitresults_shapes.hdf5 (must contain the parameter "
        "covariance, i.e. fit run with the Hessian)",
    )
    p.add_argument(
        "--effMCFile",
        required=True,
        help="in-situ effMC pkl (insitu_effMC_<proc>.pkl.lz4); needed to map the "
        "logit coefficients to SF = eps/e",
    )
    p.add_argument(
        "--prevSFFile",
        default=None,
        help="theta_central pkl of the iteration that produced this fit "
        "(None -> zeros, iteration 0)",
    )
    p.add_argument(
        "--tagAndProbeDir",
        default=None,
        help="dir with allEfficiencies_2D_<step>_<charge>.root for the "
        "external T&P SF overlay",
    )
    p.add_argument("--plotdir", required=True, help="output plot directory")
    p.add_argument(
        "--ptRange",
        type=float,
        nargs=2,
        default=[26.0, 70.0],
        help="Chebyshev pt window used in production (= --pt min max)",
    )
    p.add_argument(
        "--utSlices",
        type=float,
        nargs="+",
        default=[-10.0, 0.0, 20.0],
        help="uT values (GeV) at which to draw the HLT/Iso SF bands",
    )
    p.add_argument("--ylim", type=float, nargs=2, default=[0.7, 1.3])
    p.add_argument("--result", default=None, help="named fit result to read")
    p.add_argument(
        "--eoscp",
        action="store_true",
        help="copy to eos with xrdcp instead of the mount",
    )
    p.add_argument("--debug", action="store_true")
    args = p.parse_args()
    logging.setup_logger("plot_insitu_effSF", 4 if args.debug else 3)
    fitresult = get_fitresult(args.inputFile, result=args.result)
    if "cov" not in fitresult:
        raise ValueError(
            f"{args.inputFile} has no parameter covariance; rerun the fit with the "
            "Hessian (use the *_shapes.hdf5 output)."
        )
    make_plots(fitresult, args)


if __name__ == "__main__":
    main()
