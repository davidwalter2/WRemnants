import h5py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from uncertainty_tools import (
    all_mc_corrections,
    get_era_vals,
    get_mc_lumis,
    make_mutually_exclusive,
)

from wremnants.utilities.io_tools import input_tools
from wums.boostHistHelpers import (
    divideHists,
    multiplyHists,
)

mass_bin = 9
var_size = 0.01


matplotlib.rcParams.update({"font.size": 12})

file_in = "/work/submit/jbenke/WRemnants/scripts/histmakers/"
file_out = "/home/submit/jbenke/public_html/"
file_in_name = file_in + "mz_dilepton_liv_scetlib_dyturbo_CT18Z_N3p0LL_N2LO_Corr.hdf5"

ramses_slope = 0.0007
hfoc_slope = 0.0006


def make_plot(
    data_all,
    plotname,
    legend_all=["MC", "Data"],
    ylim=[],
    error_bars=True,
    ylabel="",
    colors=[],
    linestyles=[],
    file_out_modifier="",
):
    plt.clf()
    plt.tight_layout()
    plt.subplots_adjust(left=0.17)
    for i in range(len(legend_all)):
        if len(linestyles) != 0:
            if len(colors) != 0:
                data_all[i].plot1d(
                    yerr=error_bars, color=colors[i], linestyle=linestyles[i]
                )
        elif len(colors) != 0:
            data_all[i].plot1d(yerr=error_bars, color=colors[i])
        else:
            data_all[i].plot1d(yerr=error_bars)
    plt.xlim([0.01, 24])
    if len(ylim) != 0:
        plt.ylim(ylim)
    if len(ylabel) != 0:
        plt.ylabel(ylabel, fontsize=16)

    plt.xlabel("Sidereal Time [hr]", fontsize=16)
    # plt.title(plotname)
    plt.legend(legend_all)
    plt.savefig(file_out + file_out_modifier + plotname + ".png")


h5file = h5py.File(file_in_name, "r")
results = input_tools.load_results_h5py(h5file)

data_output = results["SingleMuon_2016PostVFP"]["output"]
lumi_output = results["SingleMuon_2016PostVFP"]["lumi_outout"]
MC_Zmumu = results["Zmumu_2016PostVFP"]["output"]

dtdt_data = data_output["time_mll"].get()
dtst_data = data_output["time_dtst"].get()
stst_data = data_output["time_stst"].get()
iso_data = data_output["time_iso"].get()

time_proj_low = data_output["time_proj"].get()

dtdt_prpg_BG, dtdt_prpg_BG_syst, dtdt_prpg_BG_stat = get_era_vals(
    MC_Zmumu, "dtdt", "BG"
)
dtst_prpg_BG, dtst_prpg_BG_syst, dtst_prpg_BG_stat = get_era_vals(
    MC_Zmumu, "dtst", "BG"
)
stst_prpg_BG, stst_prpg_BG_syst, stst_prpg_BG_stat = get_era_vals(
    MC_Zmumu, "stst", "BG"
)
iso_BG, iso_BG_syst, iso_BG_stat = get_era_vals(MC_Zmumu, "pass_iso", "BG", iso=True)


dtdt_prpg_H, dtdt_prpg_H_syst, dtdt_prpg_H_stat = get_era_vals(MC_Zmumu, "dtdt", "H")
dtst_prpg_H, dtst_prpg_H_syst, dtst_prpg_H_stat = get_era_vals(MC_Zmumu, "dtst", "H")
stst_prpg_H, stst_prpg_H_syst, stst_prpg_H_stat = get_era_vals(MC_Zmumu, "stst", "H")
iso_H, iso_H_syst, iso_H_stat = get_era_vals(MC_Zmumu, "pass_iso", "H", iso=True)


pass_gen = MC_Zmumu["pass_gen"].get()


### STABILITY
lumi_hfoc = lumi_output["lumi_hfoc"].get()
lumi_pcc = lumi_output["lumi_pcc"].get()
lumi_ramses = lumi_output["lumi_ramses"].get()

## LINEARITY
sbil_pcc = lumi_output["sbil_pcc"].get()
count_pcc = lumi_output["count_pcc"].get()


lumi_hfoc_nom = lumi_output["lumi_in_hfoc"].get()
lumi_pcc_nom = lumi_output["lumi_in_pcc"].get()
lumi_ramses_nom = lumi_output["lumi_in_ramses"].get()

lumi_scaling = lumi_output["lumi_nom"].get()
lumi_scaling_h = lumi_output["lumi_pre"].get()
lumi_scaling_bg = lumi_output["lumi_post"].get()


weightsum = results["Zmumu_2016PostVFP"]["weight_sum"]
cross_sec = results["Zmumu_2016PostVFP"]["dataset"]["xsec"]

#### A COUPLE FIXED QUANTITIES
nbins_mll = len(dtdt_data.axes["mll"])
nbins_time = len(dtst_data.axes["time"])
nbins_pt = len(dtst_data.axes["pt_probe"])
nbins_eta = len(dtdt_data.axes["eta_probe"])

hfoc_scaling = divideHists(lumi_hfoc, lumi_hfoc_nom)
hfoc_scaling = multiplyHists(hfoc_scaling, lumi_scaling)

pcc_scaling = divideHists(lumi_pcc, lumi_pcc_nom)
pcc_scaling = multiplyHists(pcc_scaling, lumi_scaling)

ramses_scaling = divideHists(lumi_ramses, lumi_ramses_nom)
ramses_scaling = multiplyHists(ramses_scaling, lumi_scaling)


time_proj_low = time_proj_low[{"mll": mass_bin}]
iso_H = iso_H[{"mll": mass_bin}]
dtdt_prpg_H = dtdt_prpg_H[{"mll": mass_bin}]
dtst_prpg_H = dtst_prpg_H[{"mll": mass_bin}]
stst_prpg_H = stst_prpg_H[{"mll": mass_bin}]

iso_BG = iso_BG[{"mll": mass_bin}]
dtdt_prpg_BG = dtdt_prpg_BG[{"mll": mass_bin}]
dtst_prpg_BG = dtst_prpg_BG[{"mll": mass_bin}]
stst_prpg_BG = stst_prpg_BG[{"mll": mass_bin}]


dtdt_data = dtdt_data[{"mll": mass_bin}]
dtst_data = dtst_data[{"mll": mass_bin}]
stst_data = stst_data[{"mll": mass_bin}]
iso_data = iso_data[{"mll": mass_bin}]

pass_gen = pass_gen[{"mll": mass_bin}]
prpg_all = [
    iso_H,
    dtdt_prpg_H,
    dtst_prpg_H,
    stst_prpg_H,
    iso_BG,
    dtdt_prpg_BG,
    dtst_prpg_BG,
    stst_prpg_BG,
]


pass_gen = all_mc_corrections(
    pass_gen,
    time_proj_low,
    lumi_scaling,
    weightsum,
    cross_sec,
)


lumi_hists = [lumi_scaling_h, lumi_scaling_bg]
iso_hfoc, dtdt_prpg_hfoc, dtst_prpg_hfoc, stst_prpg_hfoc = get_mc_lumis(
    prpg_all,
    time_proj_low,
    hfoc_scaling,
    lumi_hists,
    weightsum,
    cross_sec,
)

### lumi
hfoc_scaling = divideHists(lumi_hfoc, lumi_hfoc_nom)
pcc_scaling = divideHists(lumi_pcc, lumi_pcc_nom)
ramses_scaling = divideHists(lumi_ramses, lumi_ramses_nom)


## add the time axis to mc data
iso, dtdt_prpg, dtst_prpg, stst_prpg = get_mc_lumis(
    prpg_all,
    time_proj_low,
    lumi_scaling,
    lumi_hists,
    weightsum,
    cross_sec,
)


iso_data, dtdt_data, dtst_data, stst_data = make_mutually_exclusive(
    iso_data, dtdt_data, dtst_data, stst_data
)
#### should look into whether i need these


# pdb.set_trace()
def empty_hist_copy(hist):
    # creates am empty hist of the same shape
    copy_array = np.ones_like(hist.values())
    copy_hist = hist.copy()
    copy_hist.values()[...] = copy_array
    return copy_hist


make_plot(
    [pcc_scaling, hfoc_scaling, ramses_scaling],
    "lumi_ratios",
    ["PCC", "HFOC", "RAMSES"],
    [0.999, 1.001],
    False,
    "ratio of inst. lumi to PHYSICS",
    colors=["black", "red", "blue"],
    file_out_modifier="liv_uncert/lumi/2026-02-17/",
)

### this doesn't plot the fitted version, just the theoretical one
# make_plot([hfoc_scaling], "lumi_ratios", ["HFOC/nominal"], [0.9, 1.1],file_out_modifier="liv_uncert/lumi/2026-02-16/")
# make_plot([pcc_scaling], "lumi_ratios", ["PCC/nominal"], [0.999, 1.001], file_out_modifier="liv_uncert/lumi/2026-02-16/")
# make_plot([ramses_scaling], "lumi_ratios", ["RAMSES/nominal"], [0.9, 1.1], file_out_modifier="liv_uncert/lumi/2026-02-16/")
