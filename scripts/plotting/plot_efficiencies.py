import numpy as np
import matplotlib as mpl
import hist
import uncertainties as unc
from uncertainties import unumpy as unp

from wremnants import plot_tools
from utilities import logging, common
from utilities.io_tools import, input_tools, output_tools

import pdb

parser = common.plot_parser()
parser.add_argument('--axes', default=["nPU",], nargs="+", type=str, help='Specify x-axis to be used')
parser.add_argument('--ylim', default=None, nargs=2, type=float, help='Specify y axis range')
parser.add_argument('--rlim', default=None, nargs=2, type=float, help='Specify ratio axis range')
args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose)

colors = mpl.colormaps["tab10"]

res, meta, _ = input_tools.read_infile(args.infile)

s = hist.tag.Slicer()

selections_base = {
    "nPU": s[0:50:], # select events in mass range without underflow and overflow bins
    "mll": s[0:hist.overflow:hist.sum], # select events in mass range without underflow and overflow bins
    "acceptance_gen": s[2::hist.sum], # select in gen acceptance (1 for signal; 2 for inclusive; 4 for barrel-barrel)
}

selections = {**selections_base,
    "os": True, # select opposite charged candidates
    "acceptance_reco": s[1::hist.sum], # select in reco acceptance
    # "acceptance_gen": s[1::hist.sum], # select in gen acceptance (1 for signal; 2 for inclusive; 4 for barrel-barrel)
    }
# no opposite charge requirement for standalone muons
standalone_selections = {**selections_base,
    "acceptance_reco": s[1::hist.sum], # select in reco acceptance
    # "acceptance_gen": s[1::hist.sum], # select in gen acceptance (1 for signal; 2 for inclusive; 4 for barrel-barrel)
}
gen_selections = {**selections_base,
    "os": True, # select opposite charged candidates
    # "acceptance_gen": s[2::hist.sum], # select in gen acceptance (2 for inclusive; 4 for barrel-barrel)
}

def plot_1d(arrays, names, proc, xedges, ratio=True, xlabel=None, ylabel=None, rlabel="Ratio-1", percentage=True, percentage_ratio=True,
    xmin=None, xmax=None, ylim=None, rlim=None, outname=""
):
    logger.info(f"Make 1D plot for {names}")

    if not isinstance(arrays, list):
        arrays = [arrays]
        names = [names]

    if ylabel is None:
        ylabel = "a.u."

    if percentage:
        scale=100
        arrays = [a*scale for a in arrays]
        ylabel = f"{ylabel} [%]"
    else:
        scale=1

    if ylim is None:
        ymax = max([max(unp.nominal_values(a)) for a in arrays])
        ymin = min([min(unp.nominal_values(a)) for a in arrays])
        yrange = ymax - ymin
        ymin = ymin if ymin == 0 else ymin - yrange*0.3
        ymax = ymax + yrange*0.3
    else:
        ymax = ylim[1]
        ymin = ylim[0]

    if xmin is not None:
        xlim = (xmin, xmax)
    else:
        xlim = (xedges[0],xedges[-1])

    x_array = xedges[:-1] + (xedges[1:] - xedges[:-1])/2.

    figure_args = dict(x_array=x_array, xlabel=xlabel, ylabel=ylabel, cms_label=args.cmsDecor, 
        automatic_scale=False, width_scale=1.2, ylim=(ymin, ymax), xlim=xlim  
    )
    if ratio:
        if percentage_ratio:
            scale_ratio=100
            rlabel = f"{rlabel} [%]"
        arrays_ratio = [(divide(a,arrays[0])-1)*scale_ratio for a in arrays]
        if rlim is None:
            rmin = min([min(unp.nominal_values(a)) for a in arrays_ratio])
            rmax = max([max(unp.nominal_values(a)) for a in arrays_ratio])
        else:
            rmin, rmax = rlim
        fig, ax1, ax2 = plot_tools.figureWithRatio(**figure_args, rlabel=rlabel, rrange=(rmin, rmax))
    else:
        fig, ax1 = plot_tools.figure(**figure_args)

    ax1.plot(xlim, [scale,scale], color="black", linestyle="--")

    for i, a in enumerate(arrays):
        y = unp.nominal_values(a)
        err = unp.std_devs(a)

        ax1.stairs(y, xedges, color=colors(i), label=names[i])
        ax1.bar(x=xedges[:-1], height=2*err, bottom=y - err, width=np.diff(xedges), align='edge', linewidth=0, alpha=0.3, color=colors(i))

    plot_tools.addLegend(ax1, ncols=1, text_size=12)

    if ratio:
        ax2.plot(xlim, [0,0], color="black", linestyle="--")

        for i, a in enumerate(arrays_ratio):
            y = unp.nominal_values(a)
            err = unp.std_devs(a)

            ax2.stairs(y, xedges, color=colors(i), label=names[i])
            ax2.bar(x=xedges[:-1], height=2*err, bottom=y - err, width=np.diff(xedges), align='edge', linewidth=0, alpha=0.3, color=colors(i))

    if xlabel:
        outname += f"_{xlabel}"
    output_tools.make_plot_dir(args.outpath, args.outfolder)
    outpath = args.outpath+"/"+args.outfolder
    plot_name = f"hist_{outname}_{proc}"
    if args.postfix:
        plot_name += f"_{args.postfix}"
    plot_tools.save_pdf_and_png(outpath, plot_name)
    plot_tools.write_index_and_log(outpath, plot_name, args=args, analysis_meta_info=meta[0])


def divide(uarr1, uarr2, cutoff=1e-5):
    cutoff_mask = unp.nominal_values(uarr2) > cutoff
    out = uarr1.copy()
    out[(np.abs(unp.nominal_values(uarr1)) < cutoff) & (np.abs(unp.nominal_values(uarr2)) < cutoff)] = unc.ufloat(1, 0)
    return np.divide(uarr1, uarr2, out=out, where=cutoff_mask)

def get_array(histos, axis, name, selection=None, base_name="pair", postfx="mass", flow=False, uncertainties=True):
    histname = f"{base_name}_{name}_{axis}_{postfx}"
    if selection is not None:
        histo = histos[histname].get()[selection].project(axis)
    else:
        histo = histos[histname].get().project(axis)

    val = histo.values(flow=flow)
    std = histo.variances(flow=flow)**0.5*uncertainties

    edges = histo.axes[0].edges

    return unp.uarray(val, std), edges

for proc_name in ["ZmumuPostVFP"]:
    logger.info(f'Make plots for {proc_name}')

    result_base = res[proc_name]['output']

    for axis in args.axes:
        uarr_HLT_HLT, edges = get_array(result_base, axis, "HLT_HLT", selections)
        uarr_HLT_ID, edges = get_array(result_base, axis, "HLT_ID", selections)
        uarr_HLT_Glo, edges = get_array(result_base, axis, "HLT_Glo", selections)
        uarr_ID_ID, edges = get_array(result_base, axis, "ID_ID", selections)
        uarr_ID_Glo, edges = get_array(result_base, axis, "ID_Glo", selections)
        uarr_Glo_Glo, edges = get_array(result_base, axis, "Glo_Glo", selections)

        # don't consider uncertainties on gen events (difficult correlation/double counting)
        uarr_Gen, edges = get_array(result_base, axis, "Gen", gen_selections, uncertainties=False) 

        # tnp like efficiencies
        effHLT = divide(2*uarr_HLT_HLT, 2*uarr_HLT_HLT+uarr_HLT_ID)
        effID = divide(2*uarr_HLT_HLT+uarr_HLT_ID, 2*uarr_HLT_HLT+uarr_HLT_ID+uarr_HLT_Glo)

        # di muon efficiencies
        eff2HLT = divide(uarr_HLT_HLT, uarr_HLT_HLT+uarr_HLT_ID+uarr_ID_ID)
        eff2ID = divide(uarr_HLT_HLT+uarr_HLT_ID+uarr_ID_ID, uarr_HLT_HLT+uarr_HLT_ID+uarr_ID_ID+uarr_HLT_Glo+uarr_ID_Glo+uarr_Glo_Glo)
        eff2Glo = divide(uarr_HLT_HLT+uarr_HLT_ID+uarr_HLT_Glo+uarr_ID_ID+uarr_ID_Glo+uarr_Glo_Glo, uarr_Gen)

        plot_1d([eff2HLT, effHLT**2], [r"$\epsilon_{2\mathrm{HLT}}$", r"$\epsilon_{HLT}^2$"], proc_name, edges, 
            xlabel=axis, rlabel="C-1", ylim=(60, 120), rlim=(-0.5, 1), outname="dimuon_HLT")
        plot_1d([eff2ID, effID**2], [r"$\epsilon_{2\mathrm{ID}}$", r"$\epsilon_{ID}^2$"], proc_name, edges, 
            xlabel=axis, rlabel="C-1", ylim=(80, 110), rlim=(-0.1, 0.1), outname="dimuon_ID")

        # tracking efficiency
        uarr_HLT_StaFail, edges = get_array(result_base, axis, "HLT_StaFail", standalone_selections)
        uarr_HLT_StaPass, edges = get_array(result_base, axis, "HLT_StaPass", standalone_selections)

        effTrk = divide(uarr_HLT_StaPass, uarr_HLT_StaPass + uarr_HLT_StaFail)

        plot_1d([effTrk], [r"$\epsilon_{\mathrm{Trk}}$"], proc_name, edges, 
            xlabel=axis, ratio=False, ylim=(98, 101), outname="Trk")

        # standalone efficiency
        uarr_HLT_TrkFail, edges = get_array(result_base, axis, "HLT_TrkFail", selections)
        uarr_HLT_TrkPass, edges = get_array(result_base, axis, "HLT_TrkPass", selections)

        effSta = divide(uarr_HLT_TrkPass, uarr_HLT_TrkPass + uarr_HLT_TrkFail)

        plot_1d([effSta], [r"$\epsilon_{\mathrm{Sta}}$"], proc_name, edges, 
            xlabel=axis, ratio=False, ylim=(90, 105), outname="Sta")

        # global muon efficiency
        effGlo = effSta*effTrk

        # uarr_HLT_Gen, edges = get_array(result_base, axis, "HLT_Gen", gen_selections)
        # uarr_HLT_GenAnti, edges = get_array(result_base, axis, "HLT_GenAnti", gen_selections)
        # effGloTrue = divide(uarr_HLT_HLT+uarr_HLT_ID+uarr_HLT_Glo, uarr_HLT_Gen+uarr_HLT_GenAnti)

        # plot_1d([effGloTrue, effGlo], [r"$\epsilon_{\mathrm{Glo, true}}$", r"$\epsilon_{\mathrm{Glo}}$"], proc_name, edges, 
        #     xlabel=axis, ylim=(80, 110), outname="Glo")

        plot_1d([eff2Glo, effGlo**2], [r"$\epsilon_{2\mathrm{Glo}}$", r"$\epsilon_{\mathrm{Glo}}^2$"], proc_name, edges, 
            xlabel=axis, ylim=(80, 110), rlim=(-0.5, 1), outname="dimuon_Glo")

        # correlation coefficients
        cHLT = divide(effHLT**2, eff2HLT)
        cID = divide(effID**2, eff2ID)
        cGlo = divide(effGlo**2, eff2Glo)

        # closure
        # without correlation coefficients
        uarr_reco = divide((2*uarr_HLT_HLT+uarr_HLT_ID)**2, (4*uarr_HLT_HLT*(effID*effGlo)**2))
        plot_1d([uarr_Gen, uarr_reco], [r"$\mathrm{N}^\mathrm{Z}_{Gen}$", r"$\mathrm{N}^\mathrm{Z}_{Reco}$"], proc_name, edges, 
            xlabel=axis, rlabel="Closure-1", ylim=args.ylim, rlim=(-1, 1), outname="Closure", percentage=False, percentage_ratio=True)

        # with correlation coefficients
        uarr_reco_corr = uarr_reco * cHLT * cID * cGlo
        plot_1d([uarr_Gen, uarr_reco_corr], [r"$\mathrm{N}^\mathrm{Z}_{Gen}$", r"$\mathrm{N}^\mathrm{Z}_{Reco}$"], proc_name, edges, 
            xlabel=axis, rlabel="Closure-1", ylim=args.ylim, rlim=(-0.01, 0.01), outname="ClosureCorr", percentage=False, percentage_ratio=True)