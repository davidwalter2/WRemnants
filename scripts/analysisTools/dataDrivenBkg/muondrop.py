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
    stack = False
    labels = "predicted"

    if not isinstance(hp, list):
        hp = [hp]

    if len(hp) > 1:
        color = "black"
        stack = True
        colors = ["red", "blue"]
        labels = ["anti veto", "muon drop"]
    elif "veto" in name:
        color = "red"
        colors = "red"
    elif "drop" in name:
        color = "blue"
        colors = "blue"

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
            color=color,
            label="Transfer factor",
            stack=False,
            ax=ax1,
            flow="none",
            zorder=4,
        )

        outfile = name
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
        rlabel="Pred. / Obs",
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
        color=color,
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
        color=colors,
        label=labels,
        linestyle="-",
        linewidth=1,
        stack=stack,
        ax=ax1,
        yerr=True,
        binwnorm=1.0,
        flow="none",
        zorder=3,
    )

    hep.histplot(
        hh.divideHists(hh.sumHists(hp), ho),
        histtype="errorbar",
        color=color,
        ax=ax2,
        yerr=True,
        flow="none",
        zorder=3,
    )
    hep.histplot(
        hh.divideHists(ho, ho),
        histtype="step",
        color="black",
        ax=ax2,
        yerr=False,
        flow="none",
        zorder=3,
    )

    outfile = name
    if args.postfix:
        outfile += f"_{args.postfix}"

    plot_tools.addLegend(
        ax1,
        1,
    )

    plot_tools.fix_axes(
        ax1, ratio_axes, fig, yscale=args.yscale, logy=args.logy, noSci=args.noSciy
    )

    info = (
        {l: hp[i].sum() for i, l in enumerate(labels)}
        if isinstance(labels, list)
        else {labels: hp[0].sum()}
    )

    plot_tools.save_pdf_and_png(outdir, outfile)
    output_tools.write_index_and_log(
        outdir,
        outfile,
        analysis_meta_info={"data": ho.sum(), **info},
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
    lumi = 16.8 * 1000
    h = results[sample]["output"][name].get()
    scale = results[sample]["dataset"]["xsec"] * lumi / results[sample]["weight_sum"]
    h = hh.scaleHist(h, scale)
    return h


# get histograms for muon anti-veto method
for sample in [
    "ZmumuPostVFP",
    "DYJetsToMuMuMass10to50PostVFP",
    "QCDmuEnrichPt15PostVFP",
]:
    ho = loadHist(sample, "veto_met")
    hp = loadHist(sample, "pred_veto_met0")
    hm = loadHist(sample, "pred_veto_met1")

    plot(
        ho[{"charge": 1}].project("pt"),
        hp.project("pt"),
        name=f"{sample}_veto_pt_plus",
        xlabel="pt",
    )
    plot(
        ho[{"charge": 0}].project("pt"),
        hm.project("pt"),
        name=f"{sample}_veto_pt_minus",
        xlabel="pt",
    )

    plot(
        ho[{"charge": 1}].project("eta"),
        hp.project("eta"),
        name=f"{sample}_veto_eta_plus",
        xlabel="eta",
    )
    plot(
        ho[{"charge": 0}].project("eta"),
        hm.project("eta"),
        name=f"{sample}_veto_eta_minus",
        xlabel="eta",
    )

    h1om = hh.unrolledHist(
        ho[{"charge": 0}].project("pt", "eta"),
        binwnorm=1.0,
    )
    h1op = hh.unrolledHist(
        ho[{"charge": 1}].project("pt", "eta"),
        binwnorm=1.0,
    )
    h1p = hh.unrolledHist(
        hp.project("pt", "eta"),
        binwnorm=1.0,
    )
    h1m = hh.unrolledHist(
        hm.project("pt", "eta"),
        binwnorm=1.0,
    )

    plot(h1op, h1p, name=f"{sample}_veto_plus")
    plot(h1om, h1m, name=f"{sample}_veto_minus")


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

hobs = [
    hh.unrolledHist(
        h.project("pt" if idx == 0 else "pt1", "eta" if idx == 0 else "eta1"),
        binwnorm=1.0,
    )
    for idx, h in enumerate(hobs)
]

hpred = [
    hh.unrolledHist(
        h.project("pt" if idx == 0 else "pt1", "eta" if idx == 0 else "eta1"),
        binwnorm=1.0,
    )
    for idx, h in enumerate(hpred)
]

for idx in (0, 1):
    hprob = hh.unrolledHist(
        hprobs[idx].project("pt" if idx == 0 else "pt1", "eta" if idx == 0 else "eta1"),
        binwnorm=None,
    )

    plot(hobs[idx], hpred[idx], hprob, "muon_drop_" + ("plus" if idx == 0 else "minus"))


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

plot(hom, hp, name="eff_corrected_muon_drop_plus")
plot(hop, hm, name="eff_corrected_muon_drop_minus")

plot(hh.addHists(h1om, hom), [h1m, hm], name=f"sum_eff_corrected_minus")
plot(hh.addHists(h1op, hop), [h1p, hp], name=f"sum_eff_corrected_plus")


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

    axis_name = vari.split("_")[-1] if vari.startswith("new") else vari.split("_")[0]

    hobs_vetos = {}
    hpred_vetos = {}
    # hobs_drops = {}
    # hpred_drops = {}
    for sample in [
        "ZmumuPostVFP",
        "DYJetsToMuMuMass10to50PostVFP",
        "QCDmuEnrichPt15PostVFP",
    ]:
        if vari.startswith("new"):
            continue

        ho = loadHist(sample, f"veto_{vari}")
        hp = loadHist(sample, f"pred_veto_{vari}0")
        hm = loadHist(sample, f"pred_veto_{vari}1")

        hpred = hh.addHists(hp, hm).project(axis_name)

        ho = ho.project(axis_name)

        if axis_name in ["relIso"]:
            edges = hpred.axes["relIso"].edges
            new_edges = np.array(
                [*edges[:10], edges[14], edges[20], edges[40], edges[-1]]
            )

            ho = hh.rebinHist(ho, "relIso", new_edges)
            hpred = hh.rebinHist(hpred, "relIso", new_edges)

        plot(ho, hpred, name=f"{sample}_veto_{vari}", xlabel=axis_name)

        hobs_vetos[sample] = ho
        hpred_vetos[sample] = hpred

        # hobs_drops[sample] = hobs_drop
        # hpred_drops[sample] = hpred_drop

    # muon drop from high mass to low mass
    hobs_drop = loadHist("DYJetsToMuMuMass10to50PostVFP", f"drop_{vari}")

    hpred_drop0 = loadHist("ZmumuPostVFP", f"pred_drop_{vari}0")
    hpred_drop1 = loadHist("ZmumuPostVFP", f"pred_drop_{vari}1")

    hobs_drop = hobs_drop.project(axis_name)
    hpred_drop = hh.addHists(hpred_drop0, hpred_drop1).project(axis_name)

    if axis_name in ["relIso"]:
        edges = hpred_drop.axes["relIso"].edges
        new_edges = np.array([*edges[:10], edges[14], edges[20], edges[40], edges[-1]])
        hobs_drop = hh.rebinHist(hobs_drop, "relIso", new_edges)
        hpred_drop = hh.rebinHist(hpred_drop, "relIso", new_edges)

    plot(hobs_drop, hpred_drop, name=f"high2low_drop_{vari}", xlabel=axis_name)

    if not vari.startswith("new"):
        plot(
            hh.addHists(hobs_drop, hobs_vetos["DYJetsToMuMuMass10to50PostVFP"]),
            [hpred_vetos["DYJetsToMuMuMass10to50PostVFP"], hpred_drop],
            name=f"sum_{vari}",
            xlabel=axis_name,
        )
    # plot(hh.addHists(hobs_var_high, ho_high), [hpred_high, hpred_var_high], name=f"high_sum_{vari}", xlabel=vari)


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
