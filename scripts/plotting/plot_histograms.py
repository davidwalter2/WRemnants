import numpy as np
import hist
import mplhep as hep

from datasets.datagroups import Datagroups

hep.style.use(hep.style.ROOT)

from wremnants import plot_tools
from utilities import logging, common, boostHistHelpers as hh
from utilities.io_tools import, input_tools, output_tools

import pdb

parser = common.plot_parser()
parser.add_argument('--hists', default=["nPV",], nargs="+", type=str, help='Specify histogram to be used')
parser.add_argument('--ss', action="store_true",  help='Make same sign selection')
parser.add_argument('--axes', default=["mll",], nargs="+", type=str, help='Specify x-axis to be used')
parser.add_argument('--ylim', default=None, nargs=2, type=float, help='Specify y axis range')
parser.add_argument('--rlim', default=(0,2), nargs=2, type=float, help='Specify ratio axis range')
args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose)

res, meta, _ = input_tools.read_infile(args.infile)

def plot_h1d(hists, names, colors, labels, axis, ratio=True, normalize=False, 
    xlabel=None, ylabel="Events / bin", xlim=None, ylim=None, rlim=None, rlabel="Ratio", 
    flow=False, density=False, no_fill=False, outname=""
):
    logger.info(f"Make 1D plot for {names} with axis {axis}")

    if not isinstance(hists, list):
        hists = [hists]
        names = [names]
    
    h1ds = [h.project(axis) for h in hists]

    stack = [h for h, n in zip(hists,names) if n != "Data"]
    if "Data" in names:
        data_hist = [h for h, n in zip(hists,names) if n == "Data"][0]
        has_data=True
    else:
        has_data=False

    if flow:
        xedges = plot_tools.extendEdgesByFlow(h1ds[0])
    else:
        xedges = h1ds[0].axes.edges[0]

    if normalize:
        h1ds = [h/np.sum(h.values(flow=flow)) for h in h1ds]
    if density:
        for i, h1d in enumerate(h1ds):
            binwidths = xedges[1:]-xedges[:-1]
            hh.scaleHist(h1d, 1./binwidths, createNew=False)

    if ylim is None:
        ymax = max(np.array([h.values(flow=flow) for h in stack]).sum(axis=0))
        if has_data:
            ymax = max(ymax, max(data_hist.values(flow=flow)))
        ymin = 0
        yrange = ymax - ymin
        ymin = ymin if ymin == 0 else ymin - yrange*0.3
        ymax = ymax + yrange*0.3
        ylim = (ymin, ymax)
    if xlim is None:
        xlim = (xedges[0],xedges[-1])

    if ylabel is None:
        ylabel = "a.u."


    figure_args = dict(href=h1ds[0], xlabel=xlabel, ylabel=ylabel, cms_label=args.cmsDecor, 
        automatic_scale=False, width_scale=1.2, ylim=(ymin, ymax), xlim=xlim  
    )

    if ratio and has_data:
        fig, ax1, ax2 = plot_tools.figureWithRatio(**figure_args, rlabel=rlabel, rrange=rlim)
    else:
        fig, ax1 = plot_tools.figure(**figure_args)

    hep.histplot(
        stack,
        histtype="fill" if not no_fill else "step",
        color=colors,
        label=labels,
        stack=True,
        ax=ax1,
        zorder=1,
        flow='none',
    )

    if has_data:
        hep.histplot(
            data_hist, 
            yerr=True,
            histtype="errorbar",
            color=histInfo["Data"].color,
            label=histInfo["Data"].label,
            ax=ax1,
            flow='none',
        )

        if ratio:
            ax2.plot(xlim, [1,1], color="black", linestyle="--")

            hep.histplot(
                hh.divideHists(data_hist, sum(stack), cutoff=0.01, rel_unc=True, flow=False, by_ax_name=False),
                histtype="errorbar",
                color=histInfo["Data"].color,
                label=histInfo["Data"].label,
                yerr=True,
                linewidth=2,
                ax=ax2,
                flow='none',
            )

    plot_tools.addLegend(ax1, ncols=1, text_size=12)

    if xlabel:
        outname += f"_{xlabel}"
    output_tools.make_plot_dir(args.outpath, args.outfolder)
    outpath = args.outpath+"/"+args.outfolder
    plot_name = f"hist_{outname}"
    if args.postfix:
        plot_name += f"_{args.postfix}"
    plot_tools.save_pdf_and_png(outpath, plot_name)
    plot_tools.write_index_and_log(outpath, plot_name, args=args, analysis_meta_info=meta[0])

    
datagroups = Datagroups(args.infile)
# add out of acceptance Zmumu background
datagroups.setGenAxes("acceptance_gen")
datagroups.defineSignalBinsUnfolding("Zmumu", histToReadAxes="pair_Gen_nPU_mass")

histInfo = datagroups.groups

# merge gen bins 2,3,4 into Zmumu group
histInfo[f"Zmumu"].deleteMembers(histInfo[f"Zmumu"].members)
for i in (2,3,4):
    histInfo[f"Zmumu"].addMembers(histInfo[f"Zmumu_acceptance_gen{i}"].members, histInfo[f"Zmumu_acceptance_gen{i}"].memberOp)
    datagroups.deleteGroup(f"Zmumu_acceptance_gen{i}")

histInfo["Zmumu_acceptance_gen0"].color = "grey"
histInfo["Zmumu_acceptance_gen0"].label = r"Z$\to\mu\mu$ unmatched"
histInfo["Zmumu_acceptance_gen1"].color = "deepskyblue"
histInfo["Zmumu_acceptance_gen1"].label = r"Z$\to\mu\mu$ OOA"

# put unmatched to the bottom
dataset_names = [n for n in datagroups.getNames() 
    if n not in ["Zmumu_acceptance_gen0", "Wmunu"]] + ["Wmunu", "Zmumu_acceptance_gen0"]

s = hist.tag.Slicer()

for histo in args.hists:
    for name in (
        "HLT_HLT",
        "HLT_ID",
        "HLT_Glo",
        "HLT_HLT",
        "HLT_StaFail",
        "HLT_StaPass",
        "HLT_TrkFail",
        "HLT_TrkPass",
    ):  
        if "_Sta" in name:
            datagroups.setSelectOp(lambda h: h[{histo: s[::hist.sum], "acceptance_reco": s[1::hist.sum]}], processes=[d for d in dataset_names if d.startswith("Zmumu")])
            datagroups.setSelectOp(lambda h: h[{histo: s[::hist.sum]}], processes=[d for d in dataset_names if not d.startswith("Zmumu")])
        else:
            datagroups.setSelectOp(lambda h: h[{histo: s[::hist.sum], "os": not args.ss, "acceptance_reco": s[1::hist.sum]}], processes=[d for d in dataset_names if d.startswith("Zmumu")])
            datagroups.setSelectOp(lambda h: h[{histo: s[::hist.sum], "os": not args.ss}], processes=[d for d in dataset_names if not d.startswith("Zmumu")])

        logger.info(f'Make histogram plot for {histo} in category {name}')
        base_name = f"pair_{name}_{histo}_mass"

        datagroups.loadHistsForDatagroups(base_name, syst="", procsToRead=dataset_names)

        names = [n for n in dataset_names if histInfo[n].hists[base_name]][::-1]
        histos = [histInfo[n].hists[base_name] for n in names]
        colors = [histInfo[n].color for n in names if n != "Data"]
        labels = [histInfo[n].label for n in names if n != "Data"]

        for axis in args.axes:
            plot_h1d(histos, names, colors=colors, labels=labels, axis=axis, ratio=True, outname=base_name, rlim=args.rlim)