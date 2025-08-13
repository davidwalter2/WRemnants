import h5py
import hist
import mplhep as hep
import numpy as np

from utilities import parsing
from utilities.io_tools import input_tools
from utilities.styles import styles
from wums import boostHistHelpers as hh
from wums import logging, output_tools, plot_tools

parser = parsing.plot_parser()
parser.add_argument(
    "infile",
    help="Output file of the analysis stage, containing ND boost histograms",
)
parser.add_argument(
    "-n",
    "--baseName",
    type=str,
    help="Histogram name in the file (e.g., 'nominal')",
    default="nominal",
)
parser.add_argument("--logy", action="store_true", help="Enable log scale for y axis")

args = parser.parse_args()
logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

outdir = output_tools.make_plot_dir(args.outpath, args.outfolder, eoscp=args.eoscp)


h5file = h5py.File(args.infile, "r")

results = input_tools.load_results_h5py(h5file)


def plot(ho, hp, hprob=None, name="", xlabel=None):
    # if xlabel is None:
    xlabel = plot_tools.get_axis_label(styles, xlabel, xlabel)

    if hprob is not None:
        fig, ax1 = plot_tools.figure(
            hprob,
            xlabel=xlabel,
            ylabel="Transfer factor",
            ylim=args.ylim,
            logy=args.logy,
        )

        hep.histplot(
            hprob,
            yerr=hprob.variances() ** 0.5,
            histtype="errorbar",
            color="black",
            label="Transfer factor",
            stack=False,
            ax=ax1,
            flow="none",
            zorder=4,
        )

        outfile = f"muondrop_transferfactor_{name}"
        if args.postfix:
            outfile += f"_{args.postfix}"

        plot_tools.fix_axes(
            ax1, None, fig, yscale=args.yscale, logy=args.logy, noSci=args.noSciy
        )

        plot_tools.save_pdf_and_png(outdir, outfile)
        output_tools.write_index_and_log(
            outdir,
            outfile,
            analysis_meta_info={},
            args=args,
        )

    fig, ax1, ratio_axes = plot_tools.figureWithRatio(
        ho,
        xlabel=xlabel,
        ylabel="Events / bin",
        rlabel="Obs./Pred.",
        # cms_label=args.cmsDecor,
        # grid=True,
        # automatic_scale=False,
        # width_scale=1.2,
        # xlim=(0, 1),
        ylim=args.ylim,
        rrange=args.rrange,
        logy=args.logy,
    )
    ax2 = ratio_axes[-1]

    hep.histplot(
        ho,
        histtype="errorbar",
        color="black",
        label="observed",
        stack=False,
        ax=ax1,
        binwnorm=1.0,
        flow="none",
        zorder=4,
    )

    hep.histplot(
        hp,
        histtype="step",
        color="blue",
        label="predicted",
        linestyle="-",
        linewidth=1,
        stack=False,
        ax=ax1,
        yerr=True,
        binwnorm=1.0,
        flow="none",
        zorder=3,
    )

    hep.histplot(
        hh.divideHists(ho, hp),
        histtype="errorbar",
        color="black",
        ax=ax2,
        yerr=True,
        flow="none",
        zorder=3,
    )
    hep.histplot(
        hh.divideHists(hp, hp),
        histtype="step",
        color="black",
        ax=ax2,
        yerr=False,
        flow="none",
        zorder=3,
    )

    outfile = f"muondrop_{name}"
    if args.postfix:
        outfile += f"_{args.postfix}"

    plot_tools.addLegend(
        ax1,
        1,
    )

    plot_tools.fix_axes(
        ax1, ratio_axes, fig, yscale=args.yscale, logy=args.logy, noSci=args.noSciy
    )

    plot_tools.save_pdf_and_png(outdir, outfile)
    output_tools.write_index_and_log(
        outdir,
        outfile,
        analysis_meta_info={},
        args=args,
    )


# hpred = {}
# hobs = {}
# hprobs = {}
# for sample in ["ZmumuPostVFP", "DYJetsToMuMuMass10to50PostVFP"]:

# hpred0 = results[sample]["output"]["pred_drop0"].get()
# hpred1 = results[sample]["output"]["pred_drop1"].get()
# hpred = [hpred0, hpred1]


def loadHist(sample, name):
    h = results[sample]["output"][name].get()
    scale = results[sample]["dataset"]["xsec"] / results[sample]["weight_sum"]
    h = hh.scaleHist(h, scale)
    return h


# get histograms for muon drop method
def get_hists_muon_drop(name1, name2=None):
    if name2 is None:
        name2 = name1
    # low mass
    h1p = loadHist("DYJetsToMuMuMass10to50PostVFP", name1)
    h1m = loadHist("DYJetsToMuMuMass10to50PostVFP", name2)

    hp = h1p[{"pt": slice(26j, 56j), "eta": slice(-2.4j, 2.4j)}]
    hp_all = hp[
        {"eta1": slice(None, None, hist.sum), "pt1": slice(None, None, hist.sum)}
    ]
    hp_in = hp[
        {
            "eta1": slice(0, hist.overflow, hist.sum),
            "pt1": slice(0, hist.overflow, hist.sum),
        }
    ]
    # observerd positive muon inside, other ourside
    hp_out = hh.addHists(hp_all, hp_in, scale2=-1)

    hm = h1m[{"pt1": slice(26j, 56j), "eta1": slice(-2.4j, 2.4j)}]
    hm_all = hm[{"eta": slice(None, None, hist.sum), "pt": slice(None, None, hist.sum)}]
    hm_in = hm[
        {
            "eta": slice(0, hist.overflow, hist.sum),
            "pt": slice(0, hist.overflow, hist.sum),
        }
    ]
    # observerd negative muon inside, other ourside
    hm_out = hh.addHists(hm_all, hm_in, scale2=-1)

    # high mass
    h2p = loadHist("ZmumuPostVFP", name1)
    h2m = loadHist("ZmumuPostVFP", name2)

    h2p = h2p[{"pt": slice(26j, 56j), "eta": slice(-2.4j, 2.4j)}]
    h2m = h2m[{"pt1": slice(26j, 56j), "eta1": slice(-2.4j, 2.4j)}]
    h2p_in = h2p[
        {
            "eta1": slice(0, hist.overflow, hist.sum),
            "pt1": slice(0, hist.overflow, hist.sum),
        }
    ]
    h2m_in = h2m[
        {
            "eta": slice(0, hist.overflow, hist.sum),
            "pt": slice(0, hist.overflow, hist.sum),
        }
    ]

    # transfer evetns from high mass to low mass
    hp_prob = hh.divideHists(hp_out.project("eta", "pt"), h2p_in.project("eta", "pt"))
    hm_prob = hh.divideHists(
        hm_out.project("eta1", "pt1"), h2m_in.project("eta1", "pt1")
    )

    hm_out_pred = hh.multiplyHists(h2m_in, hm_prob)
    hp_out_pred = hh.multiplyHists(h2p_in, hp_prob)

    hobs = [hp_out, hm_out]
    hprobs = [hp_prob, hm_prob]
    hpred = [hp_out_pred, hm_out_pred]

    return hobs, hprobs, hpred


hobs, hprobs, hpred = get_hists_muon_drop("drop")

for idx in (0, 1):
    ho = hh.unrolledHist(
        hobs[idx].project("pt" if idx == 0 else "pt1", "eta" if idx == 0 else "eta1"),
        binwnorm=1.0,
    )
    hp = hh.unrolledHist(
        hpred[idx].project("pt" if idx == 0 else "pt1", "eta" if idx == 0 else "eta1"),
        binwnorm=1.0,
    )
    hprob = hh.unrolledHist(
        hprobs[idx].project("pt" if idx == 0 else "pt1", "eta" if idx == 0 else "eta1"),
        binwnorm=None,
    )

    plot(ho, hp, hprob, "plus" if idx == 0 else "minus")


### ---

# taking into account efficiency
hobs2 = loadHist("DYJetsToMuMuMass10to50PostVFP", f"drop_met")[{"met": hist.sum}]
hpred20 = loadHist("ZmumuPostVFP", f"pred_drop_met0")[{"met": hist.sum}]
hpred21 = loadHist("ZmumuPostVFP", f"pred_drop_met1")[{"met": hist.sum}]

hom = hh.unrolledHist(
    hobs2[{"charge": 0}].project("pt", "eta"),
    binwnorm=1.0,
)
hop = hh.unrolledHist(
    hobs2[{"charge": 1}].project("pt", "eta"),
    binwnorm=1.0,
)
hp = hh.unrolledHist(
    hpred20.project("pt", "eta"),
    binwnorm=1.0,
)
hm = hh.unrolledHist(
    hpred21.project("pt", "eta"),
    binwnorm=1.0,
)

plot(hom, hp, name="eff_corrected_plus")
plot(hop, hm, name="eff_corrected_minus")

check = False  # check for perfect closure

for vari in (
    "relIso",
    "relIso_vtxAgn04All",
    "dxy",
    "dxy_bs",
    "met",
    "mt",
    "phi",
    "new_met",
    "new_mt",
    "new_phi",
):
    logger.info(f"Now at '{vari}'")
    var = vari.split("_")[-1] if vari.startswith("new") else vari
    hobs2 = loadHist("DYJetsToMuMuMass10to50PostVFP", f"drop_{var}")

    var = (
        vari.split("_")[0]
        if any(vari.startswith(x) for x in ["relIso", "dxy"])
        else var
    )
    if check and not np.all(
        np.isclose(hobs[0].values(), hobs2[{var: hist.sum, "charge": 0}].values())
    ):
        raise RuntimeError("Should be the same for 'obs[0]'")
    if check and not np.all(
        np.isclose(hobs[1].values(), hobs2[{var: hist.sum, "charge": 1}].values())
    ):
        raise RuntimeError("Should be the same for 'obs[1]'")

    hpred20 = loadHist("ZmumuPostVFP", f"pred_drop_{vari}0")
    hpred21 = loadHist("ZmumuPostVFP", f"pred_drop_{vari}1")

    if check and not np.all(
        np.isclose(hpred[0].values(), hpred20[{var: hist.sum}].values(), atol=1e-5)
    ):
        raise RuntimeError("Should be the same for 'pred[0]'")
    if check and not np.all(
        np.isclose(hpred[1].values(), hpred21[{var: hist.sum}].values(), atol=1e-5)
    ):
        raise RuntimeError("Should be the same for 'pred[1]'")

    hobs_var = hobs2.project(var)
    hpred_var = hh.addHists(hpred20, hpred21).project(var)

    if var in ["relIso"]:
        edges = hobs_var.axes["relIso"].edges
        new_edges = np.array([*edges[:10], edges[14], edges[20], edges[40], edges[-1]])
        hobs_var = hh.rebinHist(hobs_var, "relIso", new_edges)
        hpred_var = hh.rebinHist(hpred_var, "relIso", new_edges)
    plot(hobs_var, hpred_var, name=vari, xlabel=vari)


# hobs, hprobs, hpred = get_hists_muon_drop("drop_mt0", "drop_mt1")

# hobs_mt = hh.addHists(*hobs).project("mt")
# hpred_mt = hh.addHists(*hpred).project("mt")
# plot(
#     hobs_mt,
#     hpred_mt,
#     name="mt"
# )

# hobs, hprobs, hpred = get_hists_muon_drop("drop_phi0", "drop_phi1")

# hobs_phi = hh.addHists(*hobs).project("phi")
# hpred_phi = hh.addHists(*hpred).project("phi")
# plot(
#     hobs_phi,
#     hpred_phi,
#     name="phi",
#     xlabel=r"$\Delta \phi(\mu, MET)$"
# )
