import argparse

import h5py
import numpy as np
from uncertainty_tools import (
    create_variation,
    eta_phi_systematic,
    get_corrected_mc,
    luminometer_syst,
    make_mutually_exclusive,
    make_ones_hist,
    prefiring_syst,
    remove_bins,
)

from rabbit import tensorwriter
from wremnants.utilities.io_tools import input_tools
from wums.boostHistHelpers import (
    addHists,
    divideHists,
    expand_hist_by_duplicate_axes,
    expand_hist_by_duplicate_axis,
    multiplyHists,
    scaleHist,
)

parser = argparse.ArgumentParser()
args = parser.parse_args()

slope_ramses = 0.0006
slope_hfoc = 0.0007
mass_bin = 9
var_size = 0.01

background_syst_names = [
    # "ZmumuPostVFP",
    "Top",
    "Diboson",
    "GGToLL_2016PostVFP",
    "QCDmuEnrichPt15_2016PostVFP",
    "Wplusmunu_2016PostVFP",
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
    "W",
    "QG_2L",
    "QG_Lnu",
    "Ztautau",
]

######################################################################
# DATA IMPORTS #

file_in = "/work/submit/jbenke/WRemnants/scripts/histmakers/"
file_in_name = file_in + "mz_dilepton_liv_scetlib_dyturbo_CT18Z_N3p0LL_N2LO_Corr.hdf5"
h5file = h5py.File(file_in_name, "r")
results = input_tools.load_results_h5py(h5file)

data_output = results["SingleMuon_2016PostVFP"]["output"]
lumi_output = results["SingleMuon_2016PostVFP"]["lumi_outout"]
MC_Zmumu = results["Zmumu_2016PostVFP"]["output"]

dtdt_data = data_output["time_mll"].get()
dtst_data = data_output["time_dtst"].get()
stst_data = data_output["time_stst"].get()
iso_data = data_output["time_iso"].get()
nbins_mll = len(dtdt_data.axes["mll"])

iso_data, dtdt_data, dtst_data, stst_data = make_mutually_exclusive(
    iso_data,
    dtdt_data,
    dtst_data,
    stst_data,
)


dtdt_data = remove_bins(dtdt_data)

h3_data = iso_data.project("time", "mll", "pt_probe", "eta_probe")
h2_data = dtdt_data.project("time", "mll", "pt_probe", "eta_probe")
h1_data = dtst_data.project("time", "mll", "pt_probe", "eta_probe")
h0_data = stst_data.project("time", "mll", "pt_probe", "eta_probe")

time_proj_low = data_output["time_proj"].get()

### STABILITY
lumi_hfoc = lumi_output["lumi_hfoc"].get()
lumi_pcc = lumi_output["lumi_pcc"].get()
lumi_ramses = lumi_output["lumi_ramses"].get()

## LINEARITY
sbil_pcc = lumi_output["sbil_pcc"].get()
count_pcc = lumi_output["count_pcc"].get()


### these are the instantaneous luminosities
lumi_hfoc_nom = lumi_output["lumi_physics_hfoc"].get()
lumi_pcc_nom = lumi_output["lumi_physics_pcc"].get()
lumi_ramses_nom = lumi_output["lumi_physics_ramses"].get()

lumi_scaling = lumi_output["lumi_nom"].get()
lumi_scaling_h = lumi_output["lumi_pre"].get()
lumi_scaling_bg = lumi_output["lumi_post"].get()

#### A COUPLE FIXED QUANTITIES
# nbins_mll = len(dtdt_data.axes["mll"])
nbins_time = len(dtst_data.axes["time"])
nbins_pt = len(dtst_data.axes["pt_probe"])
nbins_eta = len(dtdt_data.axes["eta_probe"])

hfoc_scaling = divideHists(lumi_hfoc, lumi_hfoc_nom)
hfoc_scaling = multiplyHists(hfoc_scaling, lumi_scaling)

# lumi scaling is presumably the nominal instantaneous
pcc_scaling = divideHists(lumi_pcc, lumi_pcc_nom)
pcc_scaling = multiplyHists(pcc_scaling, lumi_scaling)

ramses_scaling = divideHists(lumi_ramses, lumi_ramses_nom)
ramses_scaling = multiplyHists(ramses_scaling, lumi_scaling)

lumi_hists = [lumi_scaling_h, lumi_scaling_bg]

avg_sbil_pcc = scaleHist(divideHists(sbil_pcc, count_pcc), 1e9)
sbil_hfoc_fit = scaleHist(avg_sbil_pcc, slope_hfoc)

sbil_ones = make_ones_hist(sbil_hfoc_fit)
sbil_hfoc_fit = addHists(sbil_hfoc_fit, sbil_ones)
hfoc_sbil = multiplyHists(sbil_hfoc_fit, lumi_scaling)

sbil_ramses_fit = scaleHist(avg_sbil_pcc, slope_ramses)
sbil_ramses_fit = addHists(sbil_ramses_fit, sbil_ones)
ramses_sbil = multiplyHists(sbil_ramses_fit, lumi_scaling)

(
    Zmumu_mc,
    Zmumu_pass_gen,
    Zmumu_prefire,
    Zmumu_stat,
    Zmumu_weightsum,
    Zmumu_cross_sec,
    Zmumu_pcc_stability,
    Zmumu_hfoc_stability,
    Zmumu_ramses_stability,
    Zmumu_hfoc_linearity,
    Zmumu_ramses_linearity,
) = get_corrected_mc(
    results,
    "Zmumu_2016PostVFP",
    time_proj_low,
    lumi_hists,
    pcc_scaling,
    hfoc_scaling,
    ramses_scaling,
    hfoc_sbil,
    ramses_sbil,
    lumi_scaling,
    luminometers=True,
)

### need to convert this to
Zmumu_iso_prefire, Zmumu_dtdt_prefire, Zmumu_dtst_prefire, Zmumu_stst_prefire = (
    Zmumu_prefire
)
Zmumu_iso_mc, Zmumu_dtdt_mc, Zmumu_dtst_mc, Zmumu_stst_mc = Zmumu_mc
pass_gen = Zmumu_pass_gen

(
    Zmumu_iso_sbil_ramses,
    Zmumu_dtdt_sbil_ramses,
    Zmumu_dtst_sbil_ramses,
    Zmumu_stst_sbil_ramses,
) = Zmumu_ramses_linearity
Zmumu_iso_pcc, Zmumu_dtdt_pcc, Zmumu_dtst_pcc, Zmumu_stst_pcc = Zmumu_pcc_stability

Zmumu_iso_hfoc, Zmumu_dtdt_hfoc, Zmumu_dtst_hfoc, Zmumu_stst_hfoc = Zmumu_hfoc_stability
Zmumu_iso_ramses, Zmumu_dtdt_ramses, Zmumu_dtst_ramses, Zmumu_stst_ramses = (
    Zmumu_ramses_stability
)
(
    Zmumu_iso_sbil_hfoc,
    Zmumu_dtdt_sbil_hfoc,
    Zmumu_dtst_sbil_hfoc,
    Zmumu_stst_sbil_hfoc,
) = Zmumu_hfoc_linearity

(
    DYJets_mc,
    DYJets_pass_gen,
    DYJets_prefire,
    DYJets_stat,
    DYJets_weightsum,
    DYJets_cross_sec,
    DYJets_pcc_stability,
    DYJets_hfoc_stability,
    DYJets_ramses_stability,
    DYJets_hfoc_linearity,
    DYJets_ramses_linearity,
) = get_corrected_mc(
    results,
    "DYJetsToMuMuMass10to50_2016PostVFP",
    time_proj_low,
    lumi_hists,
    pcc_scaling,
    hfoc_scaling,
    ramses_scaling,
    hfoc_sbil,
    ramses_sbil,
    lumi_scaling,
    luminometers=True,
)


DYJets_iso_mc, DYJets_dtdt_mc, DYJets_dtst_mc, DYJets_stst_mc = DYJets_mc

DYJets_iso_prefire, DYJets_dtdt_prefire, DYJets_dtst_prefire, DYJets_stst_prefire = (
    DYJets_prefire
)

DYJets_iso_pcc, DYJets_dtdt_pcc, DYJets_dtst_pcc, DYJets_stst_pcc = DYJets_pcc_stability
DYJets_iso_hfoc, DYJets_dtdt_hfoc, DYJets_dtst_hfoc, DYJets_stst_hfoc = (
    DYJets_hfoc_stability
)
DYJets_iso_ramses, DYJets_dtdt_ramses, DYJets_dtst_ramses, DYJets_stst_ramses = (
    DYJets_ramses_stability
)

(
    DYJets_iso_sbil_hfoc,
    DYJets_dtdt_sbil_hfoc,
    DYJets_dtst_sbil_hfoc,
    DYJets_stst_sbil_hfoc,
) = DYJets_hfoc_linearity
(
    DYJets_iso_sbil_ramses,
    DYJets_dtdt_sbil_ramses,
    DYJets_dtst_sbil_ramses,
    DYJets_stst_sbil_ramses,
) = DYJets_ramses_linearity


iso_mc = addHists(Zmumu_iso_mc, DYJets_iso_mc)
dtdt_mc = addHists(Zmumu_dtdt_mc, DYJets_dtdt_mc)
dtst_mc = addHists(Zmumu_dtst_mc, DYJets_dtst_mc)
stst_mc = addHists(Zmumu_stst_mc, DYJets_stst_mc)


iso_sbil_ramses = addHists(Zmumu_iso_sbil_ramses, DYJets_iso_sbil_ramses)
dtdt_sbil_ramses = addHists(Zmumu_dtdt_sbil_ramses, DYJets_dtdt_sbil_ramses)
dtst_sbil_ramses = addHists(Zmumu_dtst_sbil_ramses, DYJets_dtst_sbil_ramses)
stst_sbil_ramses = addHists(Zmumu_stst_sbil_ramses, DYJets_stst_sbil_ramses)

iso_sbil_hfoc = addHists(Zmumu_iso_sbil_hfoc, DYJets_iso_sbil_hfoc)
dtdt_sbil_hfoc = addHists(Zmumu_dtdt_sbil_hfoc, DYJets_dtdt_sbil_hfoc)
dtst_sbil_hfoc = addHists(Zmumu_dtst_sbil_hfoc, DYJets_dtst_sbil_hfoc)
stst_sbil_hfoc = addHists(Zmumu_stst_sbil_hfoc, DYJets_stst_sbil_hfoc)

iso_ramses = addHists(Zmumu_iso_ramses, DYJets_iso_ramses)
dtdt_ramses = addHists(Zmumu_dtdt_ramses, DYJets_dtdt_ramses)
dtst_ramses = addHists(Zmumu_dtst_ramses, DYJets_dtst_ramses)
stst_ramses = addHists(Zmumu_stst_ramses, DYJets_stst_ramses)

iso_pcc = addHists(Zmumu_iso_pcc, DYJets_iso_pcc)
dtdt_pcc = addHists(Zmumu_dtdt_pcc, DYJets_dtdt_pcc)
dtst_pcc = addHists(Zmumu_dtst_pcc, DYJets_dtst_pcc)
stst_pcc = addHists(Zmumu_stst_pcc, DYJets_stst_pcc)

iso_hfoc = addHists(Zmumu_iso_hfoc, DYJets_iso_hfoc)
dtdt_hfoc = addHists(Zmumu_dtdt_hfoc, DYJets_dtdt_hfoc)
dtst_hfoc = addHists(Zmumu_dtst_hfoc, DYJets_dtst_hfoc)
stst_hfoc = addHists(Zmumu_stst_hfoc, DYJets_stst_hfoc)

iso_prefire = addHists(Zmumu_iso_prefire, DYJets_iso_prefire)
dtdt_prefire = addHists(Zmumu_dtdt_prefire, DYJets_dtdt_prefire)
dtst_prefire = addHists(Zmumu_dtst_prefire, DYJets_dtst_prefire)
stst_prefire = addHists(Zmumu_stst_prefire, DYJets_stst_prefire)


n_masked = pass_gen.project("time", "mll", "pt_probe", "eta_probe")


dtdt_all = addHists(iso_mc, dtdt_mc)
dtst_all = addHists(dtdt_all, dtst_mc)
stst_all = addHists(dtst_all, stst_mc)
iso_eff_var = divideHists(
    iso_mc.project("time", "mll", "pt_probe", "eta_probe"),
    dtdt_all.project("time", "mll", "pt_probe", "eta_probe"),
)  ## just want this to be one

iso_eff_var.values()[:, :, 0, :] = np.ones_like(iso_eff_var[{"pt_probe": 0}].values())


hlt_eff_var = divideHists(
    dtdt_all.project("time", "mll", "pt_probe", "eta_probe"),
    dtst_all.project("time", "mll", "pt_probe", "eta_probe"),
)

id_eff_var = divideHists(
    dtst_all.project("time", "mll", "pt_probe", "eta_probe"),
    stst_all.project("time", "mll", "pt_probe", "eta_probe"),
)
dtdt_mc = remove_bins(dtdt_mc)

h3 = iso_mc.project("time", "mll", "pt_probe", "eta_probe")
h2 = dtdt_mc.project("time", "mll", "pt_probe", "eta_probe")
h1 = dtst_mc.project("time", "mll", "pt_probe", "eta_probe")
h0 = stst_mc.project("time", "mll", "pt_probe", "eta_probe")


###################################################################
## create the tensor
writer = tensorwriter.TensorWriter()
##generator channel --> MAY BE WRONG BECAUSE THIS ISN'T MUTUTALLY EXCLUSIVE
### unrolling

# n_masked_ref = unrolledHist(n_masked)
h3_data_poi = h3_data.project("time", "mll")

h3_poi = h3.project("time", "mll")


n_masked = n_masked[{"mll": mass_bin}]
pass_gen = pass_gen[{"mll": mass_bin}]


### so I think these go in first
# writer.add_channel(n_masked.axes, "ch_masked", masked=True)  ## is this still correct?
# writer.add_process((divideHists(n_masked, lumi_scaling)), "Zmumu", "ch_masked")
### this is the low bin

h3_data_poi_high = remove_bins(h3_data_poi, "mll", mass_bin + 1)
h3_poi_high = remove_bins(h3_poi, "mll", mass_bin + 1)
h3_data_poi_low = remove_bins(h3_data_poi, "mll", mass_bin + 1, False)
h3_poi_low = remove_bins(h3_poi, "mll", mass_bin + 1, False)

writer.add_channel(h3_data_poi_low.axes, "ch_iso_poi_low")
writer.add_data(h3_data_poi_low, "ch_iso_poi_low")
writer.add_process(h3_poi_low, "Zmumu", "ch_iso_poi_low", signal=True)

writer.add_channel(h3_data_poi_high.axes, "ch_iso_poi_high")
writer.add_data(h3_data_poi_high, "ch_iso_poi_high")
writer.add_process(h3_poi_high, "Zmumu", "ch_iso_poi_high", signal=True)


h3_data_central = h3_data[{"mll": mass_bin}]
h3_mc_central = h3[{"mll": mass_bin}]
h2_data_central = h2_data[{"mll": mass_bin}]
h2_mc_central = h2[{"mll": mass_bin}]
h1_data_central = h1_data[{"mll": mass_bin}]
h1_mc_central = h1[{"mll": mass_bin}]
h0_data_central = h0_data[{"mll": mass_bin}]
h0_mc_central = h0[{"mll": mass_bin}]


h3_data_eff = h3_data
h3_mc_eff = h3
h2_data_eff = h2_data
h2_mc_eff = h2
h1_data_eff = h1_data
h1_mc_eff = h1
h0_data_eff = h0_data
h0_mc_eff = h0


print("begin expanding")

iso_eff_proj_pt_eta = expand_hist_by_duplicate_axes(
    h3_mc_eff,
    ["pt_probe", "eta_probe"],
    ["pt_tag", "eta_tag"],
)
dtdt_eff_proj_pt_eta = expand_hist_by_duplicate_axes(
    h2_mc_eff,
    ["pt_probe", "eta_probe"],
    ["pt_tag", "eta_tag"],
)
dtst_eff_proj_pt_eta = expand_hist_by_duplicate_axes(
    h1_mc_eff,
    ["pt_probe", "eta_probe"],
    ["pt_tag", "eta_tag"],
)
stst_eff_proj_pt_eta = expand_hist_by_duplicate_axes(
    h0_mc_eff,
    ["pt_probe", "eta_probe"],
    ["pt_tag", "eta_tag"],
)


iso_eff_proj = expand_hist_by_duplicate_axis(
    iso_eff_proj_pt_eta,
    "time",
    "gen_time",
)
dtdt_eff_proj = expand_hist_by_duplicate_axis(
    dtdt_eff_proj_pt_eta,
    "time",
    "gen_time",
)
dtst_eff_proj = expand_hist_by_duplicate_axis(
    dtst_eff_proj_pt_eta,
    "time",
    "gen_time",
)
stst_eff_proj = expand_hist_by_duplicate_axis(
    stst_eff_proj_pt_eta,
    "time",
    "gen_time",
)

hlt_eff_var = expand_hist_by_duplicate_axes(
    hlt_eff_var,
    ["time", "pt_probe", "eta_probe"],
    ["gen_time", "pt_tag", "eta_tag"],
)
id_eff_var = expand_hist_by_duplicate_axes(
    id_eff_var,
    ["time", "pt_probe", "eta_probe"],
    ["gen_time", "pt_tag", "eta_tag"],
)
iso_eff_var = expand_hist_by_duplicate_axes(
    iso_eff_var,
    ["time", "pt_probe", "eta_probe"],
    ["gen_time", "pt_tag", "eta_tag"],
)


id_eff_var_dtdt = remove_bins(id_eff_var)
iso_eff_var_dtdt = remove_bins(iso_eff_var)
hlt_eff_var_dtdt = remove_bins(hlt_eff_var)

writer.add_channel(h3_data_central.axes, "ch_iso_eff")
writer.add_data(h3_data_central, "ch_iso_eff")
writer.add_process(h3_mc_central, "Zmumu", "ch_iso_eff", signal=True)

writer.add_channel(h2_data_central.axes, "ch_dtdt_eff")
writer.add_data(h2_data_central, "ch_dtdt_eff")
writer.add_process(h2_mc_central, "Zmumu", "ch_dtdt_eff", signal=True)

writer.add_channel(h1_data_central.axes, "ch_dtst_eff")
writer.add_data(h1_data_central, "ch_dtst_eff")
writer.add_process(h1_mc_central, "Zmumu", "ch_dtst_eff", signal=True)

writer.add_channel(h0_data_central.axes, "ch_stst_eff")
writer.add_data(h0_data_central, "ch_stst_eff")
writer.add_process(h0_mc_central, "Zmumu", "ch_stst_eff", signal=True)


pass_gen = expand_hist_by_duplicate_axis(pass_gen, "time", "gen_time")
nbins_dtdt = nbins_pt * nbins_eta * ((nbins_pt - 1) + nbins_eta + nbins_time)
nbins_h1 = nbins_pt * nbins_eta * (nbins_pt + nbins_eta + nbins_time)
### so at this point i have already selected the mass bin, need to iterate over pt, eta, time
############### EFFICIENCY LOOP ##################
for i in range(nbins_pt):  # pt
    # for i in range(0, 2):  # pt
    print(f"pt bin: {i}")
    # for j in range(4, 5):  # eta
    for j in range(nbins_eta):  # eta

        print(f"eta_bin: {j}")

        if i > 0:
            v2 = dtdt_eff_proj_pt_eta[{"pt_tag": i - 1, "eta_tag": j}]
            var2 = addHists(v2 * var_size, h2_mc_eff)

            writer.add_systematic(
                var2[{"mll": mass_bin}],
                f"n_pt{i}_eta{j}",
                "Zmumu",
                "ch_dtdt_eff",
                constrained=False,
                groups=["nz"],
            )

        v3 = iso_eff_proj_pt_eta[{"pt_tag": i, "eta_tag": j}]
        var3 = addHists(v3 * var_size, h3_mc_eff)

        v1 = dtst_eff_proj_pt_eta[{"pt_tag": i, "eta_tag": j}]
        var1 = addHists(v1 * var_size, h1_mc_eff)

        v0 = stst_eff_proj_pt_eta[{"pt_tag": i, "eta_tag": j}]
        var0 = addHists(v0 * var_size, h0_mc_eff)

        writer.add_systematic(
            var3[{"mll": mass_bin}],
            f"n_pt{i}_eta{j}",
            "Zmumu",
            "ch_iso_eff",
            constrained=False,
            groups=["nz"],
        )
        nz = v3.project("time", "mll")

        writer.add_systematic(
            remove_bins(nz, "mll", mass_bin + 1, False),
            f"n_pt{i}_eta{j}",
            "Zmumu",
            "ch_iso_poi_low",
            constrained=False,
            groups=["nz"],
        )

        writer.add_systematic(
            remove_bins(nz, "mll", mass_bin + 1),
            f"n_pt{i}_eta{j}",
            "Zmumu",
            "ch_iso_poi_high",
            constrained=False,
            groups=["nz"],
        )

        writer.add_systematic(
            var1[{"mll": mass_bin}],
            f"n_pt{i}_eta{j}",
            "Zmumu",
            "ch_dtst_eff",
            constrained=False,
            groups=["nz"],
        )

        writer.add_systematic(
            var0[{"mll": mass_bin}],
            f"n_pt{i}_eta{j}",
            "Zmumu",
            "ch_stst_eff",
            constrained=False,
            groups=["nz"],
        )

        for k in range(nbins_time):  #  time
            # for k in range(22, 23):  #  time
            if i > 0:  ## we only have 1 bin beneath 25 GeV

                ##### HLT #####
                hlt_var_tag_h3 = create_variation(
                    hlt_eff_var, iso_mc, i, j, k, nbins_h1
                )
                hlt_probe_h3 = create_variation(
                    hlt_eff_var, iso_mc, i, j, k, 0, "probe"
                )
                hlt_probe_h3 = multiplyHists(
                    scaleHist(hlt_probe_h3, var_size), h3_mc_eff
                )
                hlt_var_total_h3 = addHists(hlt_var_tag_h3, hlt_probe_h3)

                writer.add_systematic(
                    addHists(
                        hlt_var_total_h3[{"mll": mass_bin}],
                        h3_mc_eff[{"mll": mass_bin}],
                    ),
                    f"hlt_pt{i}_eta{j}_time{k}",
                    "Zmumu",
                    "ch_iso_eff",
                    constrained=False,
                    groups=["eff_trig"],
                )
                ### may need to broadcast to mll instead. tbd.

                poi_hlt_proj = addHists(
                    hlt_var_total_h3.project("time", "mll"),
                    h3_mc_eff.project("time", "mll"),
                )
                writer.add_systematic(
                    remove_bins(poi_hlt_proj, "mll", mass_bin + 1),
                    f"hlt_pt{i}_eta{j}_time{k}",
                    "Zmumu",
                    "ch_iso_poi_high",
                    constrained=False,
                    groups=["eff_trig"],
                )

                writer.add_systematic(
                    remove_bins(poi_hlt_proj, "mll", mass_bin + 1, False),
                    f"hlt_pt{i}_eta{j}_time{k}",
                    "Zmumu",
                    "ch_iso_poi_low",
                    constrained=False,
                    groups=["eff_trig"],
                )

                hlt_var_tag_h2 = create_variation(
                    hlt_eff_var_dtdt, dtdt_mc, i, j, k, nbins_dtdt, h2=True
                )
                hlt_probe_h2 = create_variation(
                    hlt_eff_var_dtdt, dtdt_mc, i, j, k, 0, "probe", h2=True
                )
                hlt_probe_h2 = multiplyHists(
                    scaleHist(hlt_probe_h2, var_size), h2_mc_eff
                )
                hlt_var_total_h2 = scaleHist(addHists(hlt_var_tag_h2, hlt_probe_h2), 2)

                writer.add_systematic(
                    addHists(
                        hlt_var_total_h2[{"mll": mass_bin}],
                        h2_mc_eff[{"mll": mass_bin}],
                    ),
                    f"hlt_pt{i}_eta{j}_time{k}",
                    "Zmumu",
                    "ch_dtdt_eff",
                    constrained=False,
                    groups=["eff_trig"],
                )

                hlt_var_tag_h1 = create_variation(
                    hlt_eff_var, dtst_mc, i, j, k, nbins_h1
                )
                hlt_probe_h1 = create_variation(
                    hlt_eff_var, dtst_mc, i, j, k, 0, "probe"
                )
                hlt_probe_h1 = multiplyHists(
                    scaleHist(hlt_probe_h1, var_size), h1_mc_eff
                )
                hlt_var_total_h1 = scaleHist(
                    addHists(hlt_var_tag_h1, -1 * hlt_probe_h1), 2
                )

                writer.add_systematic(
                    addHists(
                        hlt_var_total_h1[{"mll": mass_bin}],
                        h1_mc_eff[{"mll": mass_bin}],
                    ),
                    f"hlt_pt{i}_eta{j}_time{k}",
                    "Zmumu",
                    "ch_dtst_eff",
                    constrained=False,
                    groups=["eff_trig"],
                )

                hlt_var_tag_h0 = create_variation(
                    hlt_eff_var, stst_mc, i, j, k, nbins_h1
                )
                hlt_var_total_h0 = 2 * hlt_var_tag_h0

                writer.add_systematic(
                    addHists(
                        hlt_var_total_h0[{"mll": mass_bin}],
                        h0_mc_eff[{"mll": mass_bin}],
                    ),
                    f"hlt_pt{i}_eta{j}_time{k}",
                    "Zmumu",
                    "ch_stst_eff",
                    constrained=False,
                    groups=["eff_trig"],
                )

                #### ID ######

                id_var_tag_h2 = create_variation(
                    id_eff_var_dtdt, dtdt_mc, i, j, k, nbins_dtdt, h2=True
                )
                id_probe_h2 = create_variation(
                    id_eff_var_dtdt, dtdt_mc, i, j, k, 0, "probe", h2=True
                )
                id_probe_h2 = multiplyHists(scaleHist(id_probe_h2, var_size), h2_mc_eff)
                id_var_total_h2 = addHists(id_var_tag_h2, id_probe_h2)

                writer.add_systematic(
                    addHists(
                        id_var_total_h2[{"mll": mass_bin}], h2_mc_eff[{"mll": mass_bin}]
                    ),
                    f"id_pt{i}_eta{j}_time{k}",
                    "Zmumu",
                    "ch_dtdt_eff",
                    constrained=False,
                    groups=["eff_id"],
                )

                ### ISOLATION ### --> can only be valid when hlt is valid

                iso_var_tag_h2 = create_variation(
                    iso_eff_var_dtdt, dtdt_mc, i, j, k, nbins_dtdt, h2=True
                )
                iso_probe_h2 = create_variation(
                    iso_eff_var_dtdt, dtdt_mc, i, j, k, 0, "probe", h2=True
                )
                iso_probe_h2 = multiplyHists(
                    scaleHist(iso_probe_h2, var_size), h2_mc_eff
                )
                iso_var_total_h2 = scaleHist(
                    addHists(iso_var_tag_h2, -1 * iso_probe_h2), 2
                )

                writer.add_systematic(
                    addHists(
                        iso_var_total_h2[{"mll": mass_bin}],
                        h2_mc_eff[{"mll": mass_bin}],
                    ),
                    f"iso_pt{i}_eta{j}_time{k}",
                    "Zmumu",
                    "ch_dtdt_eff",
                    constrained=False,
                    groups=["eff_iso"],
                )

            iso_var_tag_h3 = create_variation(iso_eff_var, iso_mc, i, j, k, nbins_h1)
            iso_probe_h3 = create_variation(iso_eff_var, iso_mc, i, j, k, 0, "probe")
            iso_probe_h3 = multiplyHists(scaleHist(iso_probe_h3, var_size), h3_mc_eff)
            iso_var_total_h3 = addHists(iso_var_tag_h3, iso_probe_h3)

            writer.add_systematic(
                addHists(
                    iso_var_total_h3[{"mll": mass_bin}], h3_mc_eff[{"mll": mass_bin}]
                ),
                f"iso_pt{i}_eta{j}_time{k}",
                "Zmumu",
                "ch_iso_eff",
                constrained=False,
                groups=["eff_iso"],
            )
            poi_iso_proj = addHists(
                iso_var_total_h3.project("time", "mll"),
                h3_mc_eff.project("time", "mll"),
            )

            writer.add_systematic(
                remove_bins(poi_iso_proj, "mll", mass_bin + 1),
                f"iso_pt{i}_eta{j}_time{k}",
                "Zmumu",
                "ch_iso_poi_high",
                constrained=False,
                groups=["eff_iso"],
            )

            writer.add_systematic(
                remove_bins(poi_iso_proj, "mll", mass_bin + 1, False),
                f"iso_pt{i}_eta{j}_time{k}",
                "Zmumu",
                "ch_iso_poi_low",
                constrained=False,
                groups=["eff_iso"],
            )

            iso_var_tag_h1 = create_variation(iso_eff_var, dtst_mc, i, j, k, nbins_h1)
            iso_probe_h1 = create_variation(iso_eff_var, dtst_mc, i, j, k, 0, "probe")
            iso_probe_h1 = multiplyHists(scaleHist(iso_probe_h1, var_size), h1_mc_eff)
            iso_var_total_h1 = scaleHist(addHists(iso_var_tag_h1, -1 * iso_probe_h1), 2)

            writer.add_systematic(
                addHists(
                    iso_var_total_h1[{"mll": mass_bin}], h1_mc_eff[{"mll": mass_bin}]
                ),
                f"iso_pt{i}_eta{j}_time{k}",
                "Zmumu",
                "ch_dtst_eff",
                constrained=False,
                groups=["eff_iso"],
            )

            iso_var_tag_h0 = create_variation(iso_eff_var, stst_mc, i, j, k, nbins_h1)
            iso_var_total_h0 = 2 * iso_var_tag_h0
            writer.add_systematic(
                addHists(
                    iso_var_total_h0[{"mll": mass_bin}], h0_mc_eff[{"mll": mass_bin}]
                ),
                f"iso_pt{i}_eta{j}_time{k}",
                "Zmumu",
                "ch_stst_eff",
                constrained=False,
                groups=["eff_iso"],
            )

            #### FOR THE LOWEST PT BIN ####

            ##### NORMALIZATION #####
            ### need to add an extra mass axis to this

            # for masked channel --> IS NOT MUTUALLY EXCLUSIVE
            # v_masked = pass_gen[{"gen_time": k, "pt_tag": i, "eta_tag": j}]
            # var_masked = addHists(v_masked * var_size, n_masked)
            # cross_section_masked = divideHists(var_masked, lumi_scaling)

            # writer.add_systematic(
            #     cross_section_masked,
            #     f"n_pt{i}_eta{j}_time{k}",
            #     "Zmumu",
            #     "ch_masked",
            #     constrained=False,
            #     groups=["nz"],
            # )

            ###### ID #####
            id_var_tag_h3 = create_variation(id_eff_var, iso_mc, i, j, k, nbins_h1)
            id_probe_h3 = create_variation(id_eff_var, iso_mc, i, j, k, 0, "probe")
            id_probe_h3 = multiplyHists(scaleHist(id_probe_h3, var_size), h3_mc_eff)
            id_var_total_h3 = scaleHist(addHists(id_var_tag_h3, id_probe_h3), 1 / 2)

            writer.add_systematic(
                addHists(
                    id_var_total_h3[{"mll": mass_bin}], h3_mc_eff[{"mll": mass_bin}]
                ),
                f"id_pt{i}_eta{j}_time{k}",
                "Zmumu",
                "ch_iso_eff",
                constrained=False,
                groups=["eff_id"],
            )

            poi_id_proj = addHists(
                id_var_total_h3.project("time", "mll"),
                h3_mc_eff.project("time", "mll"),
            )
            writer.add_systematic(
                remove_bins(poi_id_proj, "mll", mass_bin + 1),
                f"id_pt{i}_eta{j}_time{k}",
                "Zmumu",
                "ch_iso_poi_high",
                constrained=False,
                groups=["eff_id"],
            )

            writer.add_systematic(
                remove_bins(poi_id_proj, "mll", mass_bin + 1, False),
                f"id_pt{i}_eta{j}_time{k}",
                "Zmumu",
                "ch_iso_poi_low",
                constrained=False,
                groups=["eff_id"],
            )

            id_var_tag_h1 = create_variation(id_eff_var, dtst_mc, i, j, k, nbins_h1)
            id_probe_h1 = create_variation(id_eff_var, dtst_mc, i, j, k, 0, "probe")
            id_probe_h1 = multiplyHists(scaleHist(id_probe_h1, var_size), h1_mc_eff)
            id_var_total_h1 = addHists(id_var_tag_h1, id_probe_h1)

            writer.add_systematic(
                addHists(
                    id_var_total_h1[{"mll": mass_bin}], h1_mc_eff[{"mll": mass_bin}]
                ),
                f"id_pt{i}_eta{j}_time{k}",
                "Zmumu",
                "ch_dtst_eff",
                constrained=False,
                groups=["eff_id"],
            )

            id_var_tag_h0 = create_variation(id_eff_var, stst_mc, i, j, k, nbins_h1)
            id_probe_h0 = create_variation(id_eff_var, stst_mc, i, j, k, 0, "probe")
            id_probe_h0 = multiplyHists(scaleHist(id_probe_h0, var_size), h0_mc_eff)
            id_var_total_h0 = addHists(id_var_tag_h0, -1 * id_probe_h0)

            writer.add_systematic(
                addHists(
                    id_var_total_h0[{"mll": mass_bin}], h0_mc_eff[{"mll": mass_bin}]
                ),
                f"id_pt{i}_eta{j}_time{k}",
                "Zmumu",
                "ch_stst_eff",
                constrained=False,
                groups=["eff_id"],
            )

poi_mll_high = expand_hist_by_duplicate_axis(h3_poi_high, "mll", "gen_mll")
poi_mll_low = expand_hist_by_duplicate_axis(h3_poi_low, "mll", "gen_mll")

for i in range(nbins_mll):
    if i != mass_bin:
        if i < mass_bin:
            writer.add_systematic(
                addHists(poi_mll_low[{"gen_mll": i}] * var_size, h3_poi_low),
                f"n_mll{i}",
                "Zmumu",
                "ch_iso_poi_low",
                constrained=False,
                groups=["nz"],
            )
        else:
            print("above mass bin")
            writer.add_systematic(
                addHists(
                    poi_mll_high[{"gen_mll": i - (mass_bin + 1)}] * var_size,
                    h3_poi_high,
                ),
                f"n_mll{i}",
                "Zmumu",
                "ch_iso_poi_high",
                constrained=False,
                groups=["nz"],
            )


# for i in range(len(background_syst_names)):
#     proc_name = background_proc[i]
#     if proc_name == "Zmumu fail gen":
#         fgen = True
#     else:
#         fgen = False
#     print("proc_name: %s" % proc_name)
#     background_syst(
#         writer,
#         results,
#         background_syst_names[i],
#         time_proj_low,
#         lumi_scaling,
#         [lumi_scaling_h, lumi_scaling_bg],
#         proc_name,
#         f"bkg_{proc_name}",
#         fail_gen=fgen,
#     )

num_etaphi = len(Zmumu_stat[0].project("etaPhiRegion").values())

for i in range(num_etaphi):
    eta_phi_systematic(
        writer,
        Zmumu_stat,
        time_proj_low,
        lumi_scaling,
        lumi_hists,
        Zmumu_weightsum,
        Zmumu_cross_sec,
        i,
        mass_bin,
    )


prefiring_syst(writer, iso_prefire, dtdt_prefire, dtst_prefire, stst_prefire)
luminometer_syst(writer, "pcc", iso_pcc, dtdt_pcc, dtst_pcc, stst_pcc, "stability")

## HFOC cross detector
luminometer_syst(
    writer,
    "hfoc",
    iso_sbil_hfoc,
    dtdt_sbil_hfoc,
    dtst_sbil_hfoc,
    stst_sbil_hfoc,
    "linearity",
)


luminometer_syst(
    writer,
    "hfoc",
    iso_hfoc,
    dtdt_hfoc,
    dtst_hfoc,
    stst_hfoc,
    "stability",
)

# #### RAMSES cross detector
luminometer_syst(
    writer,
    "ramses",
    iso_ramses,
    dtdt_ramses,
    dtst_ramses,
    stst_ramses,
    "stability",
)

### RAMSES linearity
luminometer_syst(
    writer,
    "ramses",
    iso_sbil_ramses,
    dtdt_sbil_ramses,
    dtst_sbil_ramses,
    stst_sbil_ramses,
    "linearity",
)

writer.write(outfolder="./", outfilename="liv")
