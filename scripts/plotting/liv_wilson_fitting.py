import argparse

import h5py
from uncertainty_tools import (
    get_era_vals,
    get_mc_lumis,
    make_mutually_exclusive,
    mc_scaling,
)

from rabbit import tensorwriter
from wremnants.utilities.io_tools import input_tools
from wums.boostHistHelpers import (
    addHists,
)

parser = argparse.ArgumentParser()
args = parser.parse_args()


file_in = "/work/submit/jbenke/WRemnants/scripts/histmakers/"
file_in_name = file_in + "mz_dilepton_liv_scetlib_dyturbo_all_bins.hdf5"

h5file = h5py.File(file_in_name, "r")
results = input_tools.load_results_h5py(h5file)


#### REAL DATA THINGS
data_output = results["SingleMuon_2016PostVFP"]["output"]
lumi_output = results["SingleMuon_2016PostVFP"]["lumi_outout"]
# MC_Zmumu = results["Zmumu_2016PostVFP"]["output"]

dtdt_data = data_output["time_mll"].get()
dtst_data = data_output["time_dtst"].get()
stst_data = data_output["time_stst"].get()
iso_data = data_output["time_iso"].get()

time_proj_low = data_output["time_proj"].get()

#### A COUPLE FIXED QUANTITIES
nbins_mll = len(dtdt_data.axes["mll"])
nbins_time = len(dtst_data.axes["time"])
nbins_pt = len(dtst_data.axes["pt_probe"])
nbins_eta = len(dtdt_data.axes["eta_probe"])


iso_data, dtdt_data, dtst_data, stst_data = make_mutually_exclusive(
    iso_data, dtdt_data, dtst_data, stst_data
)


lumi_scaling = lumi_output["lumi_nom"].get()
lumi_scaling_h = lumi_output["lumi_pre"].get()
lumi_scaling_bg = lumi_output["lumi_post"].get()

lumi_hists = [lumi_scaling_h, lumi_scaling_bg]


############


def get_corrected_mc(results, process, time_proj_low, lumi_hists):
    MC = results[process]["output"]

    iso_BG, iso_BG_syst, iso_BG_stat = get_era_vals(MC, "pass_iso", "BG", iso=True)
    dtdt_BG, dtdt_BG_syst, dtdt_BG_stat = get_era_vals(MC, "dtdt", "BG")
    dtst_BG, dtst_BG_syst, dtst_BG_stat = get_era_vals(MC, "dtst", "BG")
    stst_BG, stst_BG_syst, stst_BG_stat = get_era_vals(MC, "stst", "BG")

    iso_H, iso_H_syst, iso_H_stat = get_era_vals(MC, "pass_iso", "H", iso=True)
    dtdt_H, dtdt_H_syst, dtdt_H_stat = get_era_vals(MC, "dtdt", "H")
    dtst_H, dtst_H_syst, dtst_H_stat = get_era_vals(MC, "dtst", "H")
    stst_H, stst_H_syst, stst_H_stat = get_era_vals(MC, "stst", "H")
    syst = [
        iso_H_syst[{"downUpVar": 0}],
        dtdt_H_syst[{"downUpVar": 0}],
        dtst_H_syst[{"downUpVar": 0}],
        stst_H_syst[{"downUpVar": 0}],
        iso_BG_syst[{"downUpVar": 0}],
        dtdt_BG_syst[{"downUpVar": 0}],
        dtst_BG_syst[{"downUpVar": 0}],
        stst_BG_syst[{"downUpVar": 0}],
    ]
    stat = [
        iso_H_stat[{"downUpVar": 0}],
        dtdt_H_stat[{"downUpVar": 0}],
        dtst_H_stat[{"downUpVar": 0}],
        stst_H_stat[{"downUpVar": 0}],
        iso_BG_stat[{"downUpVar": 0}],
        dtdt_BG_stat[{"downUpVar": 0}],
        dtst_BG_stat[{"downUpVar": 0}],
        stst_BG_stat[{"downUpVar": 0}],
    ]

    prpg_all = [
        iso_H.project("mll"),
        dtdt_H.project("mll"),
        dtst_H.project("mll"),
        stst_H.project("mll"),
        iso_BG.project("mll"),
        dtdt_BG.project("mll"),
        dtst_BG.project("mll"),
        stst_BG.project("mll"),
    ]

    weightsum = results[process]["weight_sum"]
    xsec = results[process]["dataset"]["xsec"]

    pass_gen = MC["pass_gen"].get()

    iso, dtdt, dtst, stst = get_mc_lumis(
        prpg_all,
        time_proj_low.project("time", "mll"),
        lumi_scaling,
        lumi_hists,
        weightsum,
        xsec,
    )

    pass_gen = mc_scaling(
        pass_gen.project("mll"),
        time_proj_low,
        lumi_scaling,
        weightsum,
        xsec,
    )

    corrected_mc = [iso, dtdt, dtst, stst]
    ### will need to return all the syst and stat as well but can worry about that later

    return corrected_mc, pass_gen, syst, stat


Zmumu_mc, Zmumu_pass_gen, syst, stat = get_corrected_mc(
    results, "Zmumu_2016PostVFP", time_proj_low, lumi_hists
)

DYJets_mc, DYJets_pass_gen, syst, stat = get_corrected_mc(
    results, "DYJetsToMuMuMass10to50_2016PostVFP", time_proj_low, lumi_hists
)

n_masked = addHists(DYJets_pass_gen, Zmumu_pass_gen).project("time", "mll")

Zmumu_iso, Zmumu_dtdt, Zmumu_dtst, Zmumu_stst = Zmumu_mc

DYJets_iso, DYJets_dtdt, DYJets_dtst, DYJets_stst = DYJets_mc

iso = addHists(Zmumu_iso, DYJets_iso)
dtdt = addHists(Zmumu_dtdt, DYJets_dtdt)
dtst = addHists(Zmumu_dtst, DYJets_dtst)
stst = addHists(Zmumu_stst, DYJets_stst)
###################################################################
iso_data = iso_data.project("time", "mll")

# pdb.set_trace()
# iso_unrolled = unrolledHist(iso)
# iso_data_unrolled = unrolledHist(iso_data)

iso_unrolled = iso[{"mll": 9}]
iso_data_unrolled = iso_data[{"mll": 9}]


writer = tensorwriter.TensorWriter()
writer.add_channel(iso_unrolled.axes, "ch_iso")
writer.add_data(iso_data_unrolled, "ch_iso")
writer.add_process(iso_unrolled, "Zmumu", "ch_iso", signal=True)


writer.write(outfolder="./", outfilename="wilson")
