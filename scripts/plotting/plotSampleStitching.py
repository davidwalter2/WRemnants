"""
Boson mass of Z samples generated in adjacent mass ranges, to check that they join smoothly at the
sample boundaries. Reads the stitching histograms of w_z_gen_dists.py --stitchingHists: the mass
(preFSR or LHE) vs. the LHE boson pT, which allows a common phase space with samples that have a
generator cut on the Born pT. The lower panel shows the sum of the samples relative to a smooth fit
across the boundary. With --pdfs, also the PDF variations relative to the nominal prediction.
"""

import hist
import matplotlib.pyplot as plt
import numpy as np

from wremnants.utilities import parsing, theory_utils
from wremnants.utilities.io_tools import input_tools
from wums import logging, output_tools, plot_tools

parser = parsing.plot_parser()
parser.add_argument(
    "infiles",
    nargs="+",
    help="Output files of w_z_gen_dists.py --stitchingHists containing the processes",
)
parser.add_argument(
    "--procs", nargs="+", required=True, help="Processes, by increasing mass"
)
parser.add_argument("--labels", nargs="+", help="Legend labels of the processes")
parser.add_argument(
    "--hist",
    default="stitching_massVgen",
    choices=["stitching_massVgen", "stitching_massVlhe"],
    help="Histogram (mass definition)",
)
parser.add_argument(
    "--ptVlheMin",
    type=float,
    default=0.0,
    help="Lower cut on the LHE boson pT, must be a bin edge",
)
parser.add_argument(
    "--boundary", type=float, required=True, help="Mass at which the samples meet"
)
parser.add_argument(
    "--fitRange",
    type=float,
    nargs=2,
    help="Mass range of the smooth fit of the sum (default: the plot range)",
)
parser.add_argument(
    "--fitOrder",
    type=int,
    default=3,
    help="Order of the polynomial in log(m) fitted to log(dsigma/dm)",
)
parser.add_argument(
    "--selection",
    nargs="*",
    default=[],
    help="Range of further axes of the histograms (e.g. of histmakers with a reco selection), "
    "as axis=low:high (high may be inf); axes without a range are summed over",
)
parser.add_argument(
    "--selectionLabel", help="Description of the selection, added to the legend title"
)
parser.add_argument(
    "--lowPtFraction",
    type=float,
    help="Also plot the fraction of events with LHE boson pT below this value (a bin edge) "
    "vs. the mass for each sample",
)
parser.add_argument(
    "--maxRelStatUnc",
    type=float,
    default=0.05,
    help="In the variation plots, skip bins with a larger relative statistical "
    "uncertainty of the nominal prediction",
)
parser.add_argument(
    "--pdfs",
    nargs="*",
    default=[],
    choices=theory_utils.pdfMap.keys(),
    help="Also plot these PDF sets: central value and Hessian uncertainty relative to the nominal",
)
args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
labels = args.labels if args.labels else args.procs
if len(labels) != len(args.procs):
    raise ValueError("Need one label per process")


def find_file(proc):
    for f in args.infiles:
        if proc in input_tools.read_keys(f):
            return f
    raise ValueError(f"Process {proc} not found in {args.infiles}")


variation_axes = ["pdfVar", "alphasVar"]


def reduce_axes(h):
    """apply --selection to the axes beyond mass, LHE boson pT and variations, sum the rest"""
    ranges = dict(sel.split("=") for sel in args.selection)
    for name in h.axes.name[1:]:
        if name == "ptVlhe" or name in variation_axes:
            continue
        if name in ranges:
            low, high = (float(x) for x in ranges[name].split(":"))
            stop = hist.loc(high) if np.isfinite(high) else None
            h = h[{name: slice(hist.loc(low), stop, hist.sum)}]
        else:
            h = h[{name: hist.sum}]
    return h


def read(proc, name):
    return reduce_axes(input_tools.read_and_scale(find_file(proc), proc, name))


def mass_spectrum(h):
    """values and variances along the mass for LHE boson pT above the cut (incl. overflow)"""
    ax_pt = h.axes["ptVlhe"]
    edges = np.append(ax_pt.edges, np.inf)
    if not np.any(np.isclose(edges, args.ptVlheMin)):
        raise ValueError(f"--ptVlheMin {args.ptVlheMin} is not a bin edge of {edges}")
    ipt = int(np.argmin(np.abs(edges - args.ptVlheMin)))
    other = [n for n in h.axes.name if n not in ["ptVlhe"] and not n.startswith("mass")]
    # bring to (mass, ptVlhe, [variation]): the regular mass bins, the pT bins above the cut
    # including the overflow, and the regular bins of the variation axis (no flow bins)
    h = h.project(*[h.axes.name[0], "ptVlhe", *other])
    sel = [slice(1, -1), slice(ipt, None)]
    for ax in h.axes[2:]:
        start = 1 if ax.traits.underflow else 0
        sel.append(slice(start, start + ax.size))
    sel = tuple(sel)
    vals = h.values(flow=True)[sel].sum(axis=1)
    if h.storage_type == hist.storage.Weight:
        var = h.variances(flow=True)[sel].sum(axis=1)
    else:
        var = None
    return vals, var


hists = {}
for proc in args.procs:
    hists[proc] = read(proc, args.hist)

edges = hists[args.procs[0]].axes[0].edges
centers = 0.5 * (edges[1:] + edges[:-1])
widths = np.diff(edges)
xlim = args.xlim if args.xlim else [edges[0], edges[-1]]
in_plot = (edges[:-1] >= xlim[0]) & (edges[1:] <= xlim[1])

mass_label = {
    "stitching_massVgen": r"$m_{\ell\ell}^{\mathrm{preFSR}}$ (GeV)",
    "stitching_massVlhe": r"$m_{\ell\ell}^{\mathrm{LHE}}$ (GeV)",
}[args.hist]
pt_text = (
    rf"$p_{{\mathrm{{T}}}}^{{\ell\ell,\,\mathrm{{LHE}}}} > {args.ptVlheMin:g}$ GeV"
    if args.ptVlheMin > 0
    else "inclusive"
)
if args.selectionLabel:
    pt_text = (
        f"{args.selectionLabel}, {pt_text}"
        if args.ptVlheMin > 0
        else args.selectionLabel
    )

outdir = output_tools.make_plot_dir(args.outpath, args.outfolder, eoscp=args.eoscp)
postfix = f"_{args.postfix}" if args.postfix else ""
ptname = f"ptVlhe{args.ptVlheMin:g}"

# cross sections and the smooth fit of their sum
spectra = {}
for proc in args.procs:
    vals, var = mass_spectrum(hists[proc])
    spectra[proc] = (vals / widths, np.sqrt(var) / widths)
tot = sum(spectra[p][0] for p in args.procs)
tot_err = np.sqrt(sum(spectra[p][1] ** 2 for p in args.procs))

fit_range = args.fitRange if args.fitRange else xlim
in_fit = (edges[:-1] >= fit_range[0]) & (edges[1:] <= fit_range[1]) & (tot > 0)
coeffs = np.polyfit(
    np.log(centers[in_fit]),
    np.log(tot[in_fit]),
    args.fitOrder,
    w=tot[in_fit] / tot_err[in_fit],
)
fit = np.exp(np.polyval(coeffs, np.log(centers)))
pulls = (tot[in_fit] - fit[in_fit]) / tot_err[in_fit]
chi2, ndf = np.sum(pulls**2), int(np.sum(in_fit)) - args.fitOrder - 1
logger.info(f"Fit of the sum in {fit_range}: chi2/ndf = {chi2:.1f}/{ndf}")

href = hist.Hist(
    hist.axis.Variable(edges[in_plot.nonzero()[0][0] : in_plot.nonzero()[0][-1] + 2])
)
ymax = np.max(tot[in_plot])
ymin = np.min(tot[in_plot][tot[in_plot] > 0])
fig, ax1, ratio_axes = plot_tools.figureWithRatio(
    href,
    mass_label,
    r"$d\sigma/dm_{\ell\ell}$ (pb/GeV)",
    [ymin / 3, ymax * 30],
    "Sum / fit",
    args.rrange,
    xlim=xlim,
    logy=True,
)
ax2 = ratio_axes[-1]
for proc, label, color in zip(args.procs, labels, colors):
    y, yerr = spectra[proc]
    sel = in_plot & (y > 0)
    ax1.errorbar(
        centers[sel],
        y[sel],
        yerr=yerr[sel],
        xerr=widths[sel] / 2,
        fmt="none",
        color=color,
        label=label,
        linewidth=1.5,
    )
ax1.plot(centers[in_plot], tot[in_plot], color="black", linewidth=0.8, label="Sum")
ax1.plot(
    centers[in_plot & in_fit],
    fit[in_plot & in_fit],
    color="gray",
    linestyle="--",
    linewidth=0.8,
    label=f"Fit of the sum ($\\chi^2$/ndf = {chi2:.0f}/{ndf})",
)
for ax in (ax1, ax2):
    ax.axvline(args.boundary, color="gray", linestyle=":", linewidth=1)
sel = in_plot & in_fit
ax2.errorbar(
    centers[sel],
    tot[sel] / fit[sel],
    yerr=tot_err[sel] / fit[sel],
    xerr=widths[sel] / 2,
    fmt="none",
    color="black",
)
ax2.axhline(1, color="gray", linewidth=0.8)
plot_tools.addLegend(
    ax1, ncols=1, text_size="small", loc=args.legPos, title=pt_text, reverse=False
)
plot_tools.add_cms_decor(ax1, args.cmsDecor, data=False, lumi=None, loc=args.logoPos)
outfile = f"{args.hist}_{ptname}{postfix}"
plot_tools.save_pdf_and_png(outdir, outfile)
output_tools.write_index_and_log(
    outdir,
    outfile,
    args=args,
    analysis_meta_info={
        "Fit": {"range": fit_range, "order": args.fitOrder, "chi2": chi2, "ndf": ndf}
    },
)
plt.close(fig)


def step_xy(sel, y):
    """x and y of a step line over the selected (contiguous) bins, closing the last bin"""
    return np.append(edges[:-1][sel], edges[1:][sel][-1]), np.append(y, y[-1])


def good_bins(proc):
    """bins in the plot range with a precise enough nominal prediction"""
    nom, var = mass_spectrum(hists[proc])
    rel = np.divide(np.sqrt(var), nom, out=np.full_like(nom, np.inf), where=nom > 0)
    return nom, in_plot & (nom > 0) & (rel < args.maxRelStatUnc)


def finish_variation_plot(fig, ax, outfile, title=None):
    ax.axvline(args.boundary, color="gray", linestyle=":", linewidth=1)
    ax.axhline(1, color="gray", linewidth=0.8)
    plot_tools.addLegend(
        ax,
        ncols=1,
        text_size="small",
        loc=args.legPos,
        title=title if title else pt_text,
        reverse=False,
    )
    plot_tools.add_cms_decor(ax, args.cmsDecor, data=False, lumi=None, loc=args.logoPos)
    plot_tools.save_pdf_and_png(outdir, outfile)
    output_tools.write_index_and_log(outdir, outfile, args=args)
    plt.close(fig)


# fraction of the events below a LHE boson pT value, e.g. to see whether a selection implies
# the generator cut on the Born pT of one of the samples
if args.lowPtFraction is not None:
    fig, ax = plot_tools.figure(
        href,
        mass_label,
        rf"Fraction with $p_{{\mathrm{{T}}}}^{{\ell\ell,\,\mathrm{{LHE}}}}$ < "
        f"{args.lowPtFraction:g} GeV",
        [0, 1.2],
        xlim=xlim,
    )
    for proc, label, color in zip(args.procs, labels, colors):
        h = hists[proc]
        pt_edges = np.append(h.axes["ptVlhe"].edges, np.inf)
        if not np.any(np.isclose(pt_edges, args.lowPtFraction)):
            raise ValueError(f"--lowPtFraction is not a bin edge of {pt_edges}")
        ilow = int(np.argmin(np.abs(pt_edges - args.lowPtFraction)))
        h = h.project(h.axes.name[0], "ptVlhe")
        v, var = h.values(flow=True)[1:-1], h.variances(flow=True)[1:-1]
        low, high = v[:, :ilow].sum(axis=1), v[:, ilow:].sum(axis=1)
        var_low, var_high = var[:, :ilow].sum(axis=1), var[:, ilow:].sum(axis=1)
        tot = low + high
        sel = in_plot & (tot > 0)
        frac = low[sel] / tot[sel]
        err = np.sqrt(var_low[sel] * (1 - frac) ** 2 + var_high[sel] * frac**2)
        ax.errorbar(
            centers[sel],
            frac,
            yerr=err / tot[sel],
            xerr=widths[sel] / 2,
            fmt="none",
            color=color,
            label=label,
        )
    ax.axhline(0, color="gray", linewidth=0.8)
    finish_variation_plot(
        fig,
        ax,
        f"{args.hist}_fracPtVlheBelow{args.lowPtFraction:g}{postfix}",
        title=args.selectionLabel,
    )

# PDF variations relative to the nominal prediction of each sample
for pdf in args.pdfs:
    info = theory_utils.pdfMap[pdf]
    pdf_name = info["name"]
    scale = info.get("scale", 1.0)
    fig, ax = plot_tools.figure(
        href, mass_label, f"{pdf_name[3:]} / nominal", args.rrange, xlim=xlim
    )
    for proc, label, color in zip(args.procs, labels, colors):
        nom, sel = good_bins(proc)
        h = read(proc, f"{args.hist}_{pdf_name}")
        vals, _ = mass_spectrum(h)
        names = list(h.axes["pdfVar"])
        central = vals[:, 0]
        if info["combine"] == "asymHessian":
            ups = [i for i, n in enumerate(names) if n.endswith("Up")]
            downs = [names.index(names[i][:-2] + "Down") for i in ups]
            unc = 0.5 * np.sqrt(np.sum((vals[:, ups] - vals[:, downs]) ** 2, axis=1))
        else:
            unc = np.sqrt(np.sum((vals[:, 1:] - central[:, None]) ** 2, axis=1))
        unc *= scale
        ratio, band = central[sel] / nom[sel], unc[sel] / nom[sel]
        x, y_lo = step_xy(sel, ratio - band)
        _, y_hi = step_xy(sel, ratio + band)
        ax.fill_between(x, y_lo, y_hi, step="post", alpha=0.3, color=color, label=label)
        ax.step(*step_xy(sel, ratio), where="post", color=color, linewidth=1)
    finish_variation_plot(fig, ax, f"{args.hist}_{ptname}_{pdf_name}{postfix}")

    # alpha_s variations relative to the central alpha_s of the set
    if "alphasRange" not in info:
        continue
    asr = info["alphasRange"]
    fig, ax = plot_tools.figure(
        href,
        mass_label,
        rf"{pdf_name[3:]} $\alpha_s$ $\pm$ 0.{asr} / $\alpha_s$ = 0.118",
        args.rrange,
        xlim=xlim,
    )
    for proc, label, color in zip(args.procs, labels, colors):
        h = read(proc, f"{args.hist}_{pdf_name}alphaS{asr}")
        vals, _ = mass_spectrum(h)
        _, sel = good_bins(proc)
        names_as = list(h.axes["alphasVar"])
        for ivar, linestyle in [(2, "-"), (1, "--")]:
            ax.step(
                *step_xy(sel, vals[sel, ivar] / vals[sel, 0]),
                where="post",
                color=color,
                linestyle=linestyle,
                linewidth=1.2,
                label=label if ivar == 2 else None,
            )
    as_up, as_down = (f"0.{names_as[i][3:]}" for i in (2, 1))
    finish_variation_plot(
        fig,
        ax,
        f"{args.hist}_{ptname}_{pdf_name}alphaS{postfix}",
        title=f"{pt_text}; solid: {as_up}, dashed: {as_down}",
    )
