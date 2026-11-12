import h5py
import matplotlib.pyplot as plt
import numpy as np

from wremnants.utilities import parsing
from wremnants.utilities.io_tools import input_tools
from wums import output_tools, plot_tools
from wums.boostHistHelpers import (
    addHists,
    divideHists,
)

parser = parsing.plot_parser()
parser.add_argument(
    "--procFilters",
    type=str,
    nargs="*",
    default="Zmumu",
    help="Filter to plot (default no filter, only specify if you want a subset",
)
parser.add_argument(
    "--axes",
    type=str,
    nargs="+",
    default=["pt-ptGen", "abs(eta)-absEtaGen"],
    help="Define for which axes the response matrix to be plotted",
)
parser.add_argument(
    "-n",
    "--baseName",
    type=str,
    help="Histogram base name in the file (e.g., 'nominal')",
    default="nominal",
)
parser.add_argument(
    "--histName",
    type=str,
    help="Histogram name in the file (e.g., 'nominal')",
    default="nominal",
)
parser.add_argument(
    "-c",
    "--channels",
    type=str,
    nargs="+",
    choices=["plus", "minus", "all"],
    default=["all"],
    help="Select channel to plot",
)


parser.add_argument(
    "--title",
    default="Rabbit",
    type=str,
    help="Title to be printed in upper left",
)
parser.add_argument(
    "--subtitle",
    default="",
    type=str,
    help="Subtitle to be printed after title",
)
parser.add_argument("--titlePos", type=int, default=2, help="title position")


args = parser.parse_args()

cmap = plt.get_cmap("tab10")

# logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

outdir = output_tools.make_plot_dir(args.outpath, args.outfolder, eoscp=args.eoscp)

mass_bin = 9

background_syst_names = [
    # "ZmumuPostVFP",
    "Top",
    "Diboson",
    "GGToLL_2016PostVFP",
    "QCDmuEnrichPt15_2016PostVFP",
    # "Wplusmunu_2016PostVFP",
    "Wminusmunu_2016PostVFP",
    "QGToDYQTo2L_2016PostVFP",
    "QGToWQToLNu_2016PostVFP",
    "Ztautau_2016PostVFP",
]
background_proc = [
    # "Zmumu fail gen",
    "Top",
    "Diboson",
    "GG",
    "QCD",
    # "W_plus",
    "W_minus",
    "QG_2L",
    "QG_Lnu",
    "Ztautau",
]

######################################################################
# DATA IMPORTS #

file_in = "/work/submit/jbenke/WRemnants/scripts/histmakers/"
# if not args.randTime and not args.sameSignMuon:
file_in_name = (
    file_in + "mz_dilepton_liv_scetlib_dyturbo_CT18Z_N3p0LL_N2LO_Corr_response_mat.hdf5"
)

h5file = h5py.File(file_in_name, "r")
results = input_tools.load_results_h5py(h5file)
data_output = results["SingleMuon_2016PostVFP"]["output"]
lumi_output = results["SingleMuon_2016PostVFP"]["lumi_outout"]
MC_Zmumu = results["Zmumu_2016PostVFP"]["output"]
MC_DY = results["DYJetsToMuMuMass10to50_2016PostVFP"]["output"]


h_data_upper = MC_Zmumu["mass_response"].get()
weightsum = results["Zmumu_2016PostVFP"]["weight_sum"]
cross_sec = results["Zmumu_2016PostVFP"]["dataset"]["xsec"]
h_data_upper /= weightsum
h_data_upper *= cross_sec
h_data_upper *= 1000

h_data_lower = MC_DY["mass_response"].get()
weightsum = results["DYJetsToMuMuMass10to50_2016PostVFP"]["weight_sum"]
cross_sec = results["DYJetsToMuMuMass10to50_2016PostVFP"]["dataset"]["xsec"]
h_data_upper /= weightsum
h_data_upper *= cross_sec
h_data_upper *= 1000

h_data = addHists(h_data_upper, h_data_lower)
normalization = h_data.copy()
edges = h_data.axes["recoMass"].edges
single_ax = np.array([edges[i + 1] - edges[i] for i in range(len(edges) - 1)])
single_ax = single_ax[..., np.newaxis]
normalization.values()[...] = single_ax.T * single_ax
h_data_norm = divideHists(h_data, normalization)

x = h_data.axes[0].edges
y = h_data.axes[0].edges
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_aspect("auto")

mesh = ax.pcolormesh(x, y, h_data_norm.values().T, cmap="viridis")
fig.colorbar(mesh, label="Events/Gev^2", ax=ax)

ax.set_xlabel("$m_{gen}$(GeV)")
ax.set_ylabel("$m_{reco}$ (GeV)")
plot_tools.add_decor(
    ax,
    args.title,
    args.subtitle,
    lumi=16.8,  # if args.dataName == "Data" and not args.noData else None,
    loc=args.titlePos,
    text_size=args.legSize,
)

outdir = output_tools.make_plot_dir(args.outpath, eoscp=args.eoscp)
outfile = "response_matrix"
plot_tools.save_pdf_and_png(outdir, outfile)


# stability is fraction of evens in generator bin that are observed in a reconstructed bin
# pdb.set_trace()
stability_denom = h_data.project("genMass").values()
stability_denom_exp = np.array([stability_denom for i in range(len(stability_denom))]).T
stability_denom_hist = divideHists(h_data, h_data)
stability_denom_hist.values()[...] = stability_denom_exp

stability = divideHists(h_data, stability_denom_hist)
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_aspect("auto")
z = stability.values()
mesh = ax.pcolormesh(x, y, z.T, cmap="viridis")
fig.colorbar(mesh, label="Stability", ax=ax)

ax.set_xlabel("$m_{gen}$(GeV)")
ax.set_ylabel("$m_{reco}$ (GeV)")
plot_tools.add_decor(
    ax,
    args.title,
    args.subtitle,
    lumi=16.8,  # if args.dataName == "Data" and not args.noData else None,
    loc=args.titlePos,
    text_size=args.legSize,
)

for i in range(z.shape[0]):
    half_x = (x[i + 1] - x[i]) / 2
    for j in range(z.shape[1]):
        half_y = (y[j + 1] - y[j]) / 2
        if f"{z[i, j]:.3f}" != "0.000":
            ax.text(
                x[i] + half_x,  # x position (cell center)
                y[j] + half_y,  # y position (cell center)
                f"{z[i, j]:.3f}",  # text
                ha="center",
                va="center",
                color="white",
                size=6,
            )


outdir = output_tools.make_plot_dir(args.outpath, eoscp=args.eoscp)
outfile = "stability"
plot_tools.save_pdf_and_png(outdir, outfile)


# stability is fraction of evens in generator bin that are observed in a reconstructed bin
purity_denom = h_data.project("recoMass").values()
purity_denom_exp = np.array([purity_denom for i in range(len(purity_denom))])
purity_denom_hist = divideHists(h_data, h_data)
purity_denom_hist.values()[...] = purity_denom_exp

purity = divideHists(h_data, purity_denom_hist)
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_aspect("auto")
z = purity.values()
mesh = ax.pcolormesh(x, y, z.T, cmap="viridis")
fig.colorbar(mesh, label="Purity", ax=ax)

ax.set_xlabel("$m_{gen}$(GeV)")
ax.set_ylabel("$m_{reco}$ (GeV)")
plot_tools.add_decor(
    ax,
    args.title,
    args.subtitle,
    lumi=16.8,  # if args.dataName == "Data" and not args.noData else None,
    loc=args.titlePos,
    text_size=args.legSize,
)

for i in range(z.shape[0]):
    half_x = (x[i + 1] - x[i]) / 2
    for j in range(z.shape[1]):
        half_y = (y[j + 1] - y[j]) / 2
        if f"{z[i, j]:.3f}" != "0.000":
            ax.text(
                x[i] + half_x,  # x position (cell center)
                y[j] + half_y,  # y position (cell center)
                f"{z[i, j]:.3f}",  # text
                ha="center",
                va="center",
                color="white",
                size=6,
            )


outdir = output_tools.make_plot_dir(args.outpath, eoscp=args.eoscp)
outfile = "purity"
plot_tools.save_pdf_and_png(outdir, outfile)

######################################################################3
h_data_other_bins = h_data.copy()
for i in range(0, len(h_data.values()[0])):
    h_data_other_bins.values()[i, i] = -1

fig, ax = plt.subplots(figsize=(8, 6))
ax.set_aspect("auto")

mesh = ax.pcolormesh(x, y, h_data_other_bins.values().T, cmap="viridis")

ax.set_xlabel("$m_{gen}$(GeV)")
ax.set_ylabel("$m_{reco}$ (GeV)")
plot_tools.add_decor(
    ax,
    args.title,
    args.subtitle,
    lumi=16.8,  # if args.dataName == "Data" and not args.noData else None,
    loc=args.titlePos,
    text_size=args.legSize,
)
fig.colorbar(mesh, label="Events/Gev^2", ax=ax)

outdir = output_tools.make_plot_dir(args.outpath, eoscp=args.eoscp)
outfile = "response_matrix_clipped_central_bin"
plot_tools.save_pdf_and_png(outdir, outfile)


MC_Ztautau = results["Ztautau_2016PostVFP"]["output"]

h_data = MC_Ztautau["mass_response"].get()
weightsum = results["Ztautau_2016PostVFP"]["weight_sum"]
cross_sec = results["Ztautau_2016PostVFP"]["dataset"]["xsec"]
h_data /= weightsum
h_data *= cross_sec
h_data *= 1000


normalization = h_data.copy()
edges = h_data.axes["recoMass"].edges
single_ax = np.array([edges[i + 1] - edges[i] for i in range(len(edges) - 1)])
single_ax = single_ax[..., np.newaxis]
normalization.values()[...] = single_ax.T * single_ax
h_data_norm = divideHists(h_data, normalization)

fig, ax = plt.subplots(figsize=(8, 6))
ax.set_aspect("auto")

mesh = ax.pcolormesh(x, y, h_data_norm.values().T, cmap="viridis")
fig.colorbar(mesh, label="Events/Gev^2", ax=ax)

ax.set_xlabel("$m_{gen}$(GeV)")
ax.set_ylabel("$m_{reco}$ (GeV)")

plot_tools.add_decor(
    ax,
    args.title,
    args.subtitle,
    lumi=16.8,  # if args.dataName == "Data" and not args.noData else None,
    loc=args.titlePos,
    text_size=args.legSize,
)

outdir = output_tools.make_plot_dir(args.outpath, eoscp=args.eoscp)
outfile = "z_tautau_response"
plot_tools.save_pdf_and_png(outdir, outfile)


# stability is fraction of evens in generator bin that are observed in a reconstructed bin
# pdb.set_trace()
stability_denom = h_data.project("genMass").values()
stability_denom_exp = np.array([stability_denom for i in range(len(stability_denom))]).T
stability_denom_hist = divideHists(h_data, h_data)
stability_denom_hist.values()[...] = stability_denom_exp

stability = divideHists(h_data, stability_denom_hist)
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_aspect("auto")
z = stability.values()
mesh = ax.pcolormesh(x, y, z.T, cmap="viridis")
fig.colorbar(mesh, label="Stability", ax=ax)

ax.set_xlabel("$m_{gen}$(GeV)")
ax.set_ylabel("$m_{reco}$ (GeV)")
plot_tools.add_decor(
    ax,
    args.title,
    args.subtitle,
    lumi=16.8,  # if args.dataName == "Data" and not args.noData else None,
    loc=args.titlePos,
    text_size=args.legSize,
)

for i in range(z.shape[0]):
    half_x = (x[i + 1] - x[i]) / 2
    for j in range(z.shape[1]):
        half_y = (y[j + 1] - y[j]) / 2
        if f"{z[i, j]:.3f}" != "0.000":
            ax.text(
                x[i] + half_x,  # x position (cell center)
                y[j] + half_y,  # y position (cell center)
                f"{z[i, j]:.3f}",  # text
                ha="center",
                va="center",
                color="white",
                size=6,
            )


outdir = output_tools.make_plot_dir(args.outpath, eoscp=args.eoscp)
outfile = "z_tautau_stability"
plot_tools.save_pdf_and_png(outdir, outfile)


# stability is fraction of evens in generator bin that are observed in a reconstructed bin

purity_denom = h_data.project("recoMass").values()
purity_denom_exp = np.array([purity_denom for i in range(len(purity_denom))])
purity_denom_hist = divideHists(h_data, h_data)
purity_denom_hist.values()[...] = purity_denom_exp

purity = divideHists(h_data, purity_denom_hist)
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_aspect("auto")
z = purity.values()
mesh = ax.pcolormesh(x, y, z.T, cmap="viridis")
fig.colorbar(mesh, label="Purity", ax=ax)

ax.set_xlabel("$m_{gen}$(GeV)")
ax.set_ylabel("$m_{reco}$ (GeV)")
plot_tools.add_decor(
    ax,
    args.title,
    args.subtitle,
    lumi=16.8,  # if args.dataName == "Data" and not args.noData else None,
    loc=args.titlePos,
    text_size=args.legSize,
)

for i in range(z.shape[0]):
    half_x = (x[i + 1] - x[i]) / 2
    for j in range(z.shape[1]):
        half_y = (y[j + 1] - y[j]) / 2
        if f"{z[i, j]:.3f}" != "0.000":
            ax.text(
                x[i] + half_x,  # x position (cell center)
                y[j] + half_y,  # y position (cell center)
                f"{z[i, j]:.3f}",  # text
                ha="center",
                va="center",
                color="white",
                size=6,
            )


outdir = output_tools.make_plot_dir(args.outpath, eoscp=args.eoscp)
outfile = "z_tautau_purity"
plot_tools.save_pdf_and_png(outdir, outfile)


# analysis_meta_info = None
# analysis_meta_info = {"AnalysisOutput": meta["meta_info"]}
# output_tools.write_index_and_log(
#     outdir,
#     outfile,
#     analysis_meta_info={
#         **analysis_meta_info,
#         },
#         args=args,
#     )
