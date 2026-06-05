import argparse
import os

import hist
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

from wremnants.postprocessing.datagroups.datagroups import Datagroups
from wremnants.utilities import common
from wremnants.utilities.io_tools import base_io
from wums import logging, output_tools, plot_tools

hep.style.use(hep.style.ROOT)

# Produce the MC efficiency (effMC) input for the in-situ muon efficiency
# method from the 4-category probe spectra emitted by mz_dilepton.py
# (effMCprobe_{nominal,failIso,failHLT,failID}, MC only).
#
#   nominal : 2HLT, probe passes ID & HLT & Iso
#   failIso : 2HLT, probe passes ID & HLT, fails Iso
#   failHLT : 1HLT, probe passes ID, fails HLT
#   failID  : 1HLT, probe fails ID
#
# Per-step MC efficiency (raw per-bin), computed independently for each listed
# process from its own probe spectra:
#   eff_iso     =  nominal                       / (nominal + failIso)
#   eff_trigger = (nominal + failIso)            / (nominal + failIso + failHLT)
#   eff_idip    = (nominal + failIso + failHLT)  / (nominal + failIso + failHLT
#                                                              + failID)
#
# Differential structure (matches the SMP-23-002 parameterisation variables):
#   idip   : (eta, pt, charge)        -- charge-dependent, no uT
#   trigger: (eta, pt, charge, uT)    -- charge-dependent, 2D (pt, uT)
#   iso    : (eta, pt, uT)            -- charge-inclusive, 2D (pt, uT)
#
# One output pkl is written per process: insitu_effMC_<process>[_<postfix>].pkl.lz4.
# Diagnostic plots overlay all listed processes for a side-by-side comparison.
# W processes (Wmunu / Wtaunu) live in a separate single-muon histmaker output;
# multi-input-file support will be added when those are wired in.


# Color and label palette for the diagnostic-plot overlay across processes.
PROCESS_COLORS = {
    "Zmumu": "tab:blue",
    "Ztautau": "tab:red",
    "Wmunu": "tab:green",
    "Wtaunu": "tab:purple",
    "Other": "tab:gray",
    "PhotonInduced": "tab:orange",
    "Diboson": "tab:cyan",
    "Top": "tab:brown",
}
PROCESS_LABELS = {
    "Zmumu": r"Z $\to \mu\mu$",
    "Ztautau": r"Z $\to \tau\tau$",
    "Wmunu": r"W $\to \mu\nu$",
    "Wtaunu": r"W $\to \tau\nu$",
}

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputFile", type=str, required=True)
parser.add_argument("--outpath", type=str, default="./")
parser.add_argument(
    "-p", "--postfix", type=str, help="Postfix for output file names", default=None
)
parser.add_argument("--plotdir", type=str, help="Output directory for plots")
parser.add_argument(
    "--eoscp",
    action="store_true",
    help="Copy folder to eos with xrdcp rather than using the mount",
)
parser.add_argument(
    "--process",
    type=str,
    nargs="+",
    default=["Zmumu"],
    help="Process group(s) to measure effMC for. One output pkl per process; "
    "the diagnostic plots overlay all listed processes. Each must be a "
    "datagroup present in --inputFile (e.g. Zmumu, Ztautau).",
)
parser.add_argument(
    "--tagAndProbeDir",
    type=str,
    default=None,
    help="Directory with the external tag-and-probe efficiency files "
    "(allEfficiencies_2D_<step>_<charge>.root, e.g. "
    "wremnants-data/data/muonSF/tagAndProbe/2016). When set, the per-cell "
    "eff-vs-pT plots overlay the T&P EffMC2D and EffData2D as data points "
    "with stat error bars.",
)
parser.add_argument("--debug", action="store_true", help="Print debug output")
args = parser.parse_args()

logger = logging.setup_logger("make_insitu_effMC", 4 if args.debug else 3)

categories = ["nominal", "failIso", "failHLT", "failID"]
steps = ["idip", "trigger", "iso"]

# Both the dilepton (z_dilepton) and single-muon (w_mass) histmakers emit the
# same effMCprobe_{nominal,failIso,failHLT,failID} probe spectra (eta, pt,
# charge, genUT), so the per-step recombination below is identical for both.
supported_modes = ("z_dilepton", "w_mass")
datagroups = Datagroups(args.inputFile)
if datagroups.mode not in supported_modes:
    raise ValueError(
        f"Input mode '{datagroups.mode}' not supported; expected the output of "
        f"the dilepton or single-muon histmaker (one of {supported_modes})"
    )

for cat in categories:
    datagroups.loadHistsForDatagroups(f"effMCprobe_{cat}", syst="")

groups = datagroups.groups


def _ratio(num, den):
    out = np.zeros_like(num, dtype=float)
    np.divide(num, den, out=out, where=den > 0)
    return out


def _safe_ratio_1d(num_1d, den_1d):
    out = np.zeros_like(num_1d, dtype=float)
    np.divide(num_1d, den_1d, out=out, where=den_1d > 0)
    return out


def _eff_and_var(num, den, var_num, var_den):
    """Efficiency = num/den with the binomial-like variance for the case where
    the numerator events are a subset of the denominator (so
    Cov(num, den) = Var(num)). Propagation:

        Var(eff) = (Var(num)*(1-eff)^2 + eff^2 * Var(fail)) / den^2

    with Var(fail) = Var(den) - Var(num) >= 0. For unweighted Poisson counts
    this reduces to the familiar eff*(1-eff)/den. Inputs are arbitrary-shape
    numpy arrays; outputs match. Bins with den<=0 get eff=var=0.
    """
    eff = _ratio(num, den)
    var_fail = np.maximum(var_den - var_num, 0.0)
    var_eff = np.zeros_like(eff, dtype=float)
    safe = den > 0
    var_eff[safe] = (
        var_num[safe] * (1.0 - eff[safe]) ** 2 + eff[safe] ** 2 * var_fail[safe]
    ) / den[safe] ** 2
    return eff, var_eff


def _ylim_from_arr(arr, var=None, margin=0.02):
    """Adaptive y-axis range for an efficiency array.

    ymin = max(0, min(eff - err) - margin) over positive cells   (lower bound
                                                                  clipped at 0)
    ymax = 1.0 + margin                                          (hard cap, the
                                                                  physical
                                                                  efficiency
                                                                  ceiling)

    Bins above 1 from negative MC weights or large stat error bars are clipped
    by matplotlib at ymax; the dashed unity line at y=1 is always inside the
    visible range. Fallback (no positive entries): (0, 1+margin).
    """
    ymax_cap = 1.0 + margin
    positive = arr > 0
    if not positive.any():
        return 0.0, ymax_cap
    if var is None:
        lo = arr[positive].min()
    else:
        err = np.sqrt(np.maximum(var, 0.0))
        lo = float((arr - err)[positive].min())
    return max(0.0, float(lo) - margin), ymax_cap


def _draw_unity_line(ax):
    """Dashed grey horizontal line at y=1 to mark the physical efficiency
    ceiling. No legend entry."""
    ax.axhline(1.0, linestyle="--", color="gray", linewidth=0.8, zorder=0)


def compute_effMC(process):
    """Build the per-step effMC histograms and the intermediate count arrays
    needed by the diagnostic projections, for one process group.

    Returns a dict containing:
        axes      : (eta, pt, charge, genUT) input boost-hist axes
        n         : raw counts per category (eta, pt, charge, uT)
        pass2HLT, passID, allProbes : combined sums used by the ratios
        effMC     : {step: hist.Hist}  -- the per-step boost hists for the pkl
        eff_etapt : {step: ndarray}    -- uT-integrated eff for per-cell plots
                    shapes: (eta, pt, charge) for idip/trigger; (eta, pt) for iso
    """
    h = {cat: groups[process].hists[f"effMCprobe_{cat}"] for cat in categories}

    ref_axes = list(h["nominal"].axes)
    axis_names = tuple(a.name for a in ref_axes)
    expected = ("eta", "pt", "charge", "genUT")
    if axis_names != expected:
        raise ValueError(
            f"effMCprobe[{process}] axes {axis_names} != expected {expected}; "
            "histmaker has the wrong cols/axes ordering"
        )
    for cat in categories:
        cat_names = tuple(a.name for a in h[cat].axes)
        if cat_names != axis_names:
            raise ValueError(
                f"effMCprobe_{cat} axes for {process} {cat_names} != reference {axis_names}"
            )

    eta_ax, pt_ax, charge_ax, ut_ax = ref_axes

    # raw counts (sumw) and sumw2 per category, used both for the central
    # ratios and for the binomial-like variance propagation that drives the
    # plot error bars.
    n = {cat: h[cat].values(flow=False) for cat in categories}
    var_n = {cat: h[cat].variances(flow=False) for cat in categories}
    pass2HLT = n["nominal"] + n["failIso"]
    var_pass2HLT = var_n["nominal"] + var_n["failIso"]

    # Dilepton tag-and-probe counting: a 2HLT event provides TWO valid
    # (tag, probe) pairs (either leg can be the probe; both pass ID & trigger),
    # while a 1HLT event provides a single valid probe. The histmaker already
    # fills BOTH legs of each 2HLT event into the nominal/failIso histograms
    # (each at its own kinematics), so the category sums here are the exact
    # probe counts and need no double-counting factor:
    #   - trigger: num = pass2HLT, den = pass2HLT + failHLT
    #   - idip:    num = pass2HLT + failHLT (probes passing ID),
    #              den = num + failID
    #   - iso:     num = nominal, den = nominal + failIso (both 2HLT)
    # The two 2HLT legs always have opposite charge, so they never share a bin
    # -> each bin sees at most one leg per event -> clean Poisson sumw2.

    # Per-step (num, den) sumw and sumw2 in their natural shapes:
    #   idip/trigger -> (eta, pt, charge, uT);  iso -> (eta, pt, uT)
    step_num, step_den, step_var_num, step_var_den = {}, {}, {}, {}

    step_num["trigger"] = pass2HLT
    step_den["trigger"] = pass2HLT + n["failHLT"]
    step_var_num["trigger"] = var_pass2HLT
    step_var_den["trigger"] = var_pass2HLT + var_n["failHLT"]

    passID = pass2HLT + n["failHLT"]
    var_passID = var_pass2HLT + var_n["failHLT"]
    step_num["idip"] = passID
    step_den["idip"] = passID + n["failID"]
    step_var_num["idip"] = var_passID
    step_var_den["idip"] = var_passID + var_n["failID"]

    step_num["iso"] = n["nominal"].sum(axis=2)  # charge-inclusive
    step_den["iso"] = pass2HLT.sum(axis=2)
    step_var_num["iso"] = var_n["nominal"].sum(axis=2)
    step_var_den["iso"] = var_pass2HLT.sum(axis=2)

    # idip: uT and charge dependence collapsed -> (eta, pt, charge)
    eff_idip = _ratio(step_num["idip"].sum(axis=-1), step_den["idip"].sum(axis=-1))
    h_idip = hist.Hist(
        eta_ax, pt_ax, charge_ax, name="effMC_idip", storage=hist.storage.Double()
    )
    h_idip.view(flow=False)[...] = eff_idip

    # trigger: full 4D (eta, pt, charge, uT)
    eff_trig = _ratio(step_num["trigger"], step_den["trigger"])
    h_trig = hist.Hist(
        eta_ax,
        pt_ax,
        charge_ax,
        ut_ax,
        name="effMC_trigger",
        storage=hist.storage.Double(),
    )
    h_trig.view(flow=False)[...] = eff_trig

    # iso: charge-inclusive (eta, pt, uT)
    eff_iso = _ratio(step_num["iso"], step_den["iso"])
    h_iso = hist.Hist(
        eta_ax, pt_ax, ut_ax, name="effMC_iso", storage=hist.storage.Double()
    )
    h_iso.view(flow=False)[...] = eff_iso

    # uT-integrated eff per cell with its propagated variance, used by the
    # per-cell pt curves with error bars.
    eff_etapt, var_etapt = {}, {}
    for step in steps:
        eff_etapt[step], var_etapt[step] = _eff_and_var(
            step_num[step].sum(axis=-1),
            step_den[step].sum(axis=-1),
            step_var_num[step].sum(axis=-1),
            step_var_den[step].sum(axis=-1),
        )

    return {
        "axes": ref_axes,
        "n": n,
        "var_n": var_n,
        "step_num": step_num,
        "step_den": step_den,
        "step_var_num": step_var_num,
        "step_var_den": step_var_den,
        "effMC": {"idip": h_idip, "trigger": h_trig, "iso": h_iso},
        "eff_etapt": eff_etapt,
        "var_etapt": var_etapt,
    }


# Compute effMC for every requested process
process_data = {}
for process in args.process:
    if process not in groups:
        logger.warning(
            f"Process group '{process}' not in input; skipping "
            f"(available: {list(groups.keys())})"
        )
        continue
    process_data[process] = compute_effMC(process)
    for step, h_out in process_data[process]["effMC"].items():
        e = h_out.view(flow=False)
        inrange = (e > 0) & (e < 1)
        nfilled = int(inrange.sum())
        if nfilled:
            sel = e[inrange]
            logger.info(
                f"effMC[{process}][{step}] shape={e.shape}: {nfilled} bins in (0,1), "
                f"mean={sel.mean():.4f} min={sel.min():.4f} max={sel.max():.4f}; "
                f"{int((e >= 1).sum())} bins >=1 (low stats), "
                f"{int((e <= 0).sum())} empty"
            )
        else:
            logger.warning(
                f"effMC[{process}][{step}] shape={e.shape}: no bins in (0,1) "
                f"({int((e >= 1).sum())} bins >=1, {int((e <= 0).sum())} empty) "
                f"- expected only with very low statistics"
            )

if not process_data:
    raise ValueError(
        f"None of the requested processes {args.process} were found in input"
    )

# Reference axes (assumed consistent across processes — same histmaker output).
ref_axes = process_data[next(iter(process_data))]["axes"]
eta_ax, pt_ax, charge_ax, ut_ax = ref_axes


def _step_num_den(process, step):
    """Numerator/denominator (sumw and sumw2) arrays for one (process, step),
    with the 2HLT double counting already applied (see compute_effMC). Shapes
    (eta, pt, charge, uT) for idip/trigger, (eta, pt, uT) for iso.
    Returns (num, den, var_num, var_den, axis_names, charge_dep).
    """
    pd = process_data[process]
    if step in ("idip", "trigger"):
        axis_names, charge_dep = ("eta", "pt", "charge", "genUT"), True
    elif step == "iso":
        axis_names, charge_dep = ("eta", "pt", "genUT"), False
    else:
        raise KeyError(step)
    return (
        pd["step_num"][step],
        pd["step_den"][step],
        pd["step_var_num"][step],
        pd["step_var_den"][step],
        axis_names,
        charge_dep,
    )


def _step_proj_axes(step):
    """Projection axes to draw per step (IDIP is uT-independent by design)."""
    if step == "idip":
        return ("eta", "pt")
    return ("eta", "pt", "genUT")


def _project_eff(num, den, var_num, var_den, ax_idx, q_idx=None):
    """Project num/den (sumw) and var_num/var_den (sumw2) onto the requested
    axis, optionally keeping charge separately. Returns (eff, var_eff) on the
    projected shape.

    Sums of (sumw, sumw2) are additive across summed bins, so the variance of
    the projected num/den is just the sum of the input variances. The variance
    of the resulting efficiency is then computed by _eff_and_var.
    """
    if q_idx is not None:
        sum_axes = tuple(i for i in range(num.ndim) if i not in (ax_idx, q_idx))
        num_p = num.sum(axis=sum_axes)
        den_p = den.sum(axis=sum_axes)
        var_num_p = var_num.sum(axis=sum_axes)
        var_den_p = var_den.sum(axis=sum_axes)
        if ax_idx > q_idx:
            num_p, den_p = num_p.T, den_p.T
            var_num_p, var_den_p = var_num_p.T, var_den_p.T
    else:
        sum_axes = tuple(i for i in range(num.ndim) if i != ax_idx)
        num_p = num.sum(axis=sum_axes)
        den_p = den.sum(axis=sum_axes)
        var_num_p = var_num.sum(axis=sum_axes)
        var_den_p = var_den.sum(axis=sum_axes)
    return _eff_and_var(num_p, den_p, var_num_p, var_den_p)


def _proc_color(process):
    return PROCESS_COLORS.get(process)


def _proc_label(process):
    return PROCESS_LABELS.get(process, process)


# step -> (filename charge token); idip/trigger are charge-split, iso is "both".
_TNP_FILE_CHARGE = {
    "idip": {"minus": "minus", "plus": "plus"},
    "trigger": {"minus": "minus", "plus": "plus"},
    "iso": {"both": "both"},
}
_tnp_cache = {}


def load_tnp_eff(tnp_dir, step, charge_tag):
    """Read the external tag-and-probe EffMC2D and EffData2D for one
    (step, charge) from allEfficiencies_2D_<step>_<charge>.root.

    Returns a dict {"mc": (pt_centers, pt_halfwidths, eff[48,npt], err[48,npt]),
                    "data": (pt_centers, pt_halfwidths, eff[48,npt], err[48,npt])}
    with eta on the first axis (matching our 48-bin grid) and pT on the
    second. ``pt_halfwidths`` is half the (variable) pT bin width, used as the
    horizontal error bar so each point shows its pT bin extent. Results are
    cached per (step, charge). Returns None if the file is missing.
    """
    import ROOT

    key = (tnp_dir, step, charge_tag)
    if key in _tnp_cache:
        return _tnp_cache[key]

    file_charge = _TNP_FILE_CHARGE[step][charge_tag]
    path = os.path.join(tnp_dir, f"allEfficiencies_2D_{step}_{file_charge}.root")
    if not os.path.isfile(path):
        logger.warning(f"T&P file not found, skipping overlay: {path}")
        _tnp_cache[key] = None
        return None

    f = ROOT.TFile.Open(path)
    out = {}
    for which, hname in (("mc", "EffMC2D"), ("data", "EffData2D")):
        h2 = f.Get(hname)
        if not h2:
            logger.warning(f"{hname} not found in {path}")
            continue
        n_eta = h2.GetNbinsX()
        n_pt = h2.GetNbinsY()
        pt_centers = np.array([h2.GetYaxis().GetBinCenter(j + 1) for j in range(n_pt)])
        pt_halfwidths = np.array(
            [0.5 * h2.GetYaxis().GetBinWidth(j + 1) for j in range(n_pt)]
        )
        eff = np.empty((n_eta, n_pt))
        err = np.empty((n_eta, n_pt))
        for i in range(n_eta):
            for j in range(n_pt):
                eff[i, j] = h2.GetBinContent(i + 1, j + 1)
                err[i, j] = h2.GetBinError(i + 1, j + 1)
        out[which] = (pt_centers, pt_halfwidths, eff, err)
    f.Close()
    _tnp_cache[key] = out
    return out


# Marker styling for the T&P overlay points (distinct from the process lines).
_TNP_STYLE = {
    "mc": dict(fmt="o", color="black", label="T&P MC", markersize=4, capsize=2),
    "data": dict(fmt="s", color="dimgray", label="T&P data", markersize=4, capsize=2),
}


def make_effMC_plots(plotdir):
    """Layout:

        <plotdir>/<step>/effMC_<step>_proj_<axis>[_q<tag>].{pdf,png}
            1D projection: one figure per axis × charge (charge-dep steps) or
            per axis (iso). One curve per process, coloured by process.

        <plotdir>/<step>/per_cell/effMC_<step>_pt_eta<i>[_q<tag>].{pdf,png}
            eff(pT) at a fixed eta (and charge), one curve per process.

    Per-step shared y-axis range across (cells × processes × charges) so the
    plots are directly comparable. Dashed unity line on every figure.
    """
    outdir = output_tools.make_plot_dir(*plotdir.rsplit("/", 1), eoscp=args.eoscp)
    # Nest the input file's metadata under its basename so write_logfile's
    # `{**info, **meta_info}` merge does *not* clobber this script's own
    # command/time/args/git_hash entries -- same pattern as make_theory_corr.py.
    input_meta_raw = base_io.get_metadata(args.inputFile)
    input_meta_local = (
        {os.path.basename(args.inputFile): input_meta_raw} if input_meta_raw else {}
    )

    axis_objs = {"eta": eta_ax, "pt": pt_ax, "genUT": ut_ax}
    axis_xlabels = {
        "eta": r"$\eta_\mu$",
        "pt": r"$p_T^\mu$ [GeV]",
        "genUT": r"$u_T$ [GeV]",
    }
    charge_tags = {-1: ("minus", r"$\mu^{-}$"), +1: ("plus", r"$\mu^{+}$")}

    procs_present = list(process_data.keys())
    # Reference process whose efficiency drives the shared y-axis range. Other
    # processes (e.g. Ztautau) have noisy low-stat per-bin efficiencies that
    # would otherwise blow up the range, so they are excluded from the ylim
    # computation but still drawn.
    ref_proc = "Zmumu" if "Zmumu" in procs_present else procs_present[0]

    for step in steps:
        step_dir = os.path.join(outdir, step)
        per_cell_dir = os.path.join(step_dir, "per_cell")
        os.makedirs(step_dir, exist_ok=True)
        os.makedirs(per_cell_dir, exist_ok=True)

        # Determine step layout from the first process (consistent across procs).
        _, _, _, _, axis_names, charge_dep = _step_num_den(procs_present[0], step)

        # ---- 1D projections (per axis × per charge for charge-dep steps) ----
        for axis_name in _step_proj_axes(step):
            ax_idx = axis_names.index(axis_name)
            q_idx = axis_names.index("charge") if charge_dep else None
            x_ax = axis_objs[axis_name]

            charge_iter = list(range(charge_ax.size)) if charge_dep else [None]
            for j_q in charge_iter:
                # Compute the projection (eff + stat variance) per process.
                per_proc_eff = {}
                per_proc_var = {}
                for process in procs_present:
                    num, den, var_num, var_den, _, _ = _step_num_den(process, step)
                    eff_p, var_p = _project_eff(
                        num, den, var_num, var_den, ax_idx, q_idx
                    )
                    if j_q is not None:
                        per_proc_eff[process] = eff_p[:, j_q]
                        per_proc_var[process] = var_p[:, j_q]
                    else:
                        per_proc_eff[process] = eff_p
                        per_proc_var[process] = var_p

                # Shared y-axis range from the reference process only (incl.
                # its error bars); noisy processes are excluded from the range.
                ymin, ymax = _ylim_from_arr(
                    per_proc_eff[ref_proc], per_proc_var[ref_proc]
                )

                fig, ax = plt.subplots(figsize=(8, 6))
                bin_centers = x_ax.centers
                for process in procs_present:
                    y = per_proc_eff[process]
                    err = np.sqrt(np.maximum(per_proc_var[process], 0.0))
                    color = _proc_color(process)
                    ax.step(
                        x_ax.edges,
                        np.concatenate([y[:1], y]),
                        where="pre",
                        color=color,
                        label=_proc_label(process),
                    )
                    ax.errorbar(
                        bin_centers,
                        y,
                        yerr=err,
                        fmt="none",
                        ecolor=color,
                        elinewidth=1.0,
                        capsize=0,
                    )
                ax.set_xlabel(axis_xlabels[axis_name])
                ax.set_ylabel(f"MC efficiency ({step})")
                ax.set_ylim(ymin, ymax)
                _draw_unity_line(ax)
                ax.legend(loc="best", frameon=False)

                if j_q is not None:
                    qkey = int(round(charge_ax.centers[j_q]))
                    tag, qlabel = charge_tags.get(qkey, (str(qkey), f"q={qkey}"))
                    ax.text(
                        0.04,
                        0.06,
                        qlabel,
                        transform=ax.transAxes,
                        ha="left",
                        va="bottom",
                        fontsize=18,
                    )
                    suffix = f"_q{tag}"
                else:
                    suffix = ""

                plot_tools.add_decor(
                    ax, "CMS", "Simulation", data=False, lumi=None, loc=0
                )
                plot_name = f"effMC_{step}_proj_{axis_name}{suffix}"
                plot_tools.save_pdf_and_png(step_dir, plot_name)
                plt.close(fig)
                output_tools.write_index_and_log(
                    step_dir,
                    plot_name,
                    args=args,
                    analysis_meta_info=input_meta_local,
                )

        # ---- per-(eta, [q]) eff(pT) cells, uT integrated, overlaid procs ----
        # IDIP/trigger have a charge axis; iso does not.
        eff_is_charge_dep = process_data[procs_present[0]]["eff_etapt"][step].ndim == 3

        # T&P overlay (optional): preload EffMC2D/EffData2D per charge tag.
        # charge tag for the file lookup: "minus"/"plus" per cell for the
        # charge-dep steps, "both" for iso.
        tnp_by_tag = {}
        if args.tagAndProbeDir:
            if eff_is_charge_dep:
                tags_needed = [
                    charge_tags[int(round(charge_ax.centers[jq]))][0]
                    for jq in range(charge_ax.size)
                ]
            else:
                tags_needed = ["both"]
            for tag in tags_needed:
                tnp_by_tag[tag] = load_tnp_eff(args.tagAndProbeDir, step, tag)

        # Shared y-axis range from the reference process only (across cells ×
        # charges), plus the T&P points (eff ± err); noisy processes excluded.
        ylim_vals = process_data[ref_proc]["eff_etapt"][step].ravel()
        ylim_vars = process_data[ref_proc]["var_etapt"][step].ravel()
        for tnp in tnp_by_tag.values():
            if not tnp:
                continue
            for which in ("mc", "data"):
                if which in tnp:
                    _, _, e_tnp, err_tnp = tnp[which]
                    ylim_vals = np.concatenate([ylim_vals, e_tnp.ravel()])
                    ylim_vars = np.concatenate([ylim_vars, (err_tnp**2).ravel()])
        ymin, ymax = _ylim_from_arr(ylim_vals, ylim_vars)

        pt_centers = pt_ax.centers

        for i_eta in range(eta_ax.size):
            eta_lo, eta_hi = eta_ax.edges[i_eta], eta_ax.edges[i_eta + 1]
            eta_text = f"${eta_lo:+.2f}\\,\\leq\\,\\eta\\,<\\,{eta_hi:+.2f}$"
            cell_iter = range(charge_ax.size) if eff_is_charge_dep else [None]
            for j_q in cell_iter:
                fig, ax = plt.subplots(figsize=(8, 6))
                for process in procs_present:
                    eff_arr = process_data[process]["eff_etapt"][step]
                    var_arr = process_data[process]["var_etapt"][step]
                    if j_q is not None:
                        y = eff_arr[i_eta, :, j_q]
                        v = var_arr[i_eta, :, j_q]
                    else:
                        y = eff_arr[i_eta, :]
                        v = var_arr[i_eta, :]
                    err = np.sqrt(np.maximum(v, 0.0))
                    color = _proc_color(process)
                    ax.step(
                        pt_ax.edges,
                        np.concatenate([y[:1], y]),
                        where="pre",
                        color=color,
                        label=_proc_label(process),
                    )
                    ax.errorbar(
                        pt_centers,
                        y,
                        yerr=err,
                        fmt="none",
                        ecolor=color,
                        elinewidth=1.0,
                        capsize=0,
                    )
                if j_q is not None:
                    qkey = int(round(charge_ax.centers[j_q]))
                    tag, qlabel = charge_tags.get(qkey, (str(qkey), f"q={qkey}"))
                    annot = f"{eta_text}\n{qlabel}"
                    suffix = f"_q{tag}"
                else:
                    tag = "both"
                    annot = eta_text
                    suffix = ""

                # Overlay the external T&P efficiencies for this (step, eta, q)
                tnp = tnp_by_tag.get(tag)
                if tnp:
                    for which in ("mc", "data"):
                        if which not in tnp:
                            continue
                        pt_tnp, pt_hw_tnp, e_tnp, err_tnp = tnp[which]
                        ax.errorbar(
                            pt_tnp,
                            e_tnp[i_eta],
                            yerr=err_tnp[i_eta],
                            xerr=pt_hw_tnp,
                            **_TNP_STYLE[which],
                        )
                ax.text(
                    0.04,
                    0.06,
                    annot,
                    transform=ax.transAxes,
                    ha="left",
                    va="bottom",
                    fontsize=18,
                )
                ax.set_xlabel(r"$p_T^\mu$ [GeV]")
                ax.set_ylabel(f"MC efficiency ({step})")
                ax.set_ylim(ymin, ymax)
                _draw_unity_line(ax)
                ax.legend(loc="best", frameon=False)
                plot_tools.add_decor(
                    ax, "CMS", "Simulation", data=False, lumi=None, loc=0
                )
                plot_name = f"effMC_{step}_pt_eta{i_eta:02d}{suffix}"
                plot_tools.save_pdf_and_png(per_cell_dir, plot_name)
                plt.close(fig)
                output_tools.write_index_and_log(
                    per_cell_dir,
                    plot_name,
                    args=args,
                    analysis_meta_info=input_meta_local,
                )

    logger.info(f"Wrote effMC plots to {outdir}")
    if output_tools.is_eosuser_path(plotdir) and args.eoscp:
        output_tools.copy_to_eos(outdir, plotdir)


if args.plotdir:
    make_effMC_plots(args.plotdir)


# Write one pkl per process.
filename_base = "insitu_effMC"
if args.postfix:
    filename_base += f"_{args.postfix}"
input_meta = base_io.get_metadata(args.inputFile)
for process, pd in process_data.items():
    outfile = f"{args.outpath}/{filename_base}_{process}.pkl.lz4"
    output_dict = {"effMC": pd["effMC"], "process": process, "steps": steps}
    # Carry over the input histmaker's provenance as file_meta_data, and let
    # write_lz4_pkl_output attach this script's own meta_data (command, git
    # hash, timestamp) so the output is inspectable with
    # scripts/inspect/print_command.py.
    output_tools.write_lz4_pkl_output(
        outfile,
        process,
        output_dict,
        common.base_dir,
        args,
        input_meta,
    )
    logger.info(f"Wrote in-situ effMC[{process}] to {outfile}")
