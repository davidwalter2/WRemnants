import argparse
import time

import h5py
import hist
import numpy as np
from liv_variation_helpers_rapidity import (
    background_syst,
    eta_phi_systematic,
    get_corrected_mc,
    luminometer_syst,
    make_mutually_exclusive,
    make_ones_hist,
    poi_projection,
    prefiring_syst,
    remove_bins,
    variation_values,
)

from rabbit import tensorwriter
from wremnants.utilities.io_tools import input_tools
from wums.boostHistHelpers import (
    addHists,
    divideHists,
    expand_hist_by_duplicate_axis,
    multiplyHists,
    scaleHist,
)

parser = argparse.ArgumentParser()
parser.add_argument(
    "--randTime",
    type=lambda s: s.lower() in ("true", "1", "yes"),
    default=False,
    help="use assign times by event number",
)

parser.add_argument(
    "--sameSignMuon",
    type=lambda s: s.lower() in ("true", "1", "yes"),
    default=False,
    help="whether to select muons with the same sign in events",
)
args = parser.parse_args()

slope_ramses = 0.0006
slope_hfoc = 0.0007
mass_bin = 9
var_size = 0.01

background_syst_names = [
    # "ZmumuPostVFP",
    "Top",
    # "Diboson",
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
    # "Diboson",
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
file_in_name = file_in + "mz_dilepton_liv_scetlib_dyturbo_CT18Z_N3p0LL_N2LO_Corr.hdf5"
if args.randTime:
    file_in_name = (
        file_in + "mz_dilepton_liv_scetlib_dyturbo_CT18Z_N3p0LL_N2LO_Corr_randTime.hdf5"
    )
if args.sameSignMuon:
    file_in_name = (
        file_in
        + "mz_dilepton_liv_scetlib_dyturbo_CT18Z_N3p0LL_N2LO_Corr_sameSignMuons.hdf5"
    )
h5file = h5py.File(file_in_name, "r")
results = input_tools.base_io.load_results_h5py(h5file)
data_output = results["SingleMuon_2016PostVFP"]["output"]
lumi_output = results["SingleMuon_2016PostVFP"]["lumi_outout"]
MC_Zmumu = results["Zmumu_2016PostVFP"]["output"]

dtdt_data = data_output["time_mll"].get()
dtst_data = data_output["time_dtst"].get()
stst_data = data_output["time_stst"].get()
iso_data = data_output["time_iso"].get()

iso_data, dtdt_data, dtst_data, stst_data = make_mutually_exclusive(
    iso_data,
    dtdt_data,
    dtst_data,
    stst_data,
)
dtdt_data = remove_bins(dtdt_data)

#### A COUPLE FIXED QUANTITIES
nbins_mll = len(dtdt_data.axes["mll"])
nbins_time = len(dtst_data.axes["time"])
nbins_pt = len(dtst_data.axes["pt_probe"])
nbins_eta = len(dtdt_data.axes["eta_probe"])

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

# lumi scaling is presumably the nominal instantaneous
hfoc_scaling = multiplyHists(divideHists(lumi_hfoc, lumi_hfoc_nom), lumi_scaling)
pcc_scaling = multiplyHists(divideHists(lumi_pcc, lumi_pcc_nom), lumi_scaling)
ramses_scaling = multiplyHists(divideHists(lumi_ramses, lumi_ramses_nom), lumi_scaling)

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
    "Zmumu10to50_2016PostVFP",
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


# def sum_regions(region1, region2):
#     ### need to convert this to


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


probe_axes = ("time", "mll", "pt_probe", "eta_probe")
iso_all = iso_mc.project(*probe_axes)
dtdt_all = addHists(iso_all, dtdt_mc.project(*probe_axes))
dtst_all = addHists(dtdt_all, dtst_mc.project(*probe_axes))
stst_all = addHists(dtst_all, stst_mc.project(*probe_axes))


iso_eff_var = divideHists(iso_all, dtdt_all)  ## just want this to be one

iso_eff_var.values()[:, :, 0, :] = np.ones_like(iso_eff_var[{"pt_probe": 0}].values())

hlt_eff_var = divideHists(dtdt_all, dtst_all)

id_eff_var = divideHists(dtst_all, stst_all)
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
# POI channels are (time, mll, eta_probe), with the eta_probe binning of the efficiency channels
h3_data_poi = poi_projection(h3_data)
h3_poi = poi_projection(h3)


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


print("begin preparing variation arrays")

hlt_eff_var_values = variation_values(hlt_eff_var)
id_eff_var_values = variation_values(id_eff_var)
iso_eff_var_values = variation_values(iso_eff_var)


# efficiencies depend only on (pt, eta): use the mll = mass_bin values in every mll bin,
# so the same tag and probe variations are applied to the POI channels
def _mll_independent(values):
    return np.broadcast_to(values[:, mass_bin : mass_bin + 1], values.shape)


hlt_eff_var_values = _mll_independent(hlt_eff_var_values)
id_eff_var_values = _mll_independent(id_eff_var_values)
iso_eff_var_values = _mll_independent(iso_eff_var_values)
h2_pt_offset = nbins_pt - len(h2_mc_eff.axes["pt_probe"])
h2_eta_offset = nbins_eta - len(h2_mc_eff.axes["eta_probe"])
hlt_eff_var_dtdt_values = hlt_eff_var_values[:, :, h2_pt_offset:, h2_eta_offset:]
id_eff_var_dtdt_values = id_eff_var_values[:, :, h2_pt_offset:, h2_eta_offset:]
iso_eff_var_dtdt_values = iso_eff_var_values[:, :, h2_pt_offset:, h2_eta_offset:]


def _batched_gen_axis():
    return hist.axis.StrCategory(
        [f"time{k}" for k in range(nbins_time)], name="gen_time"
    )


def _batched_hist(values, reference_hist):
    return hist.Hist(
        _batched_gen_axis(),
        *reference_hist.axes,
        data=values,
    )


def _normalization_array(reference_values, i, j, pt_offset=0, eta_offset=0):
    i -= pt_offset
    j -= eta_offset
    output = np.zeros_like(reference_values)
    output[:, :, i, j] = reference_values[:, :, i, j] * var_size
    return output


def _tag_variation_array(
    variation_values,
    reference_tag_values,
    i,
    j,
    h2=False,
    pt_offset=0,
    eta_offset=0,
):
    """Build all gen-time tag variations with NumPy.

    Varying the efficiency in bin (i, j) changes every probe bin through the
    events whose tag muon is in (i, j), so the variation is spread over all
    probe bins using the tag-resolved MC yield.
    """
    hist_i = i - pt_offset if h2 else i
    hist_j = j - eta_offset if h2 else j
    factor = variation_values[:, :, hist_i, hist_j]
    tag_in_bin = reference_tag_values[:, :, :, :, hist_i, hist_j]
    output = np.zeros(
        (factor.shape[0], *tag_in_bin.shape),
        dtype=np.result_type(factor, tag_in_bin),
    )
    diagonal = np.arange(factor.shape[0])
    output[diagonal, diagonal] = factor[:, :, None, None] * tag_in_bin * var_size
    return output


def _probe_variation_array(
    variation_values, reference_values, i, j, pt_offset=0, eta_offset=0
):
    """Build all gen-time probe variations with NumPy."""
    probe_i = i - pt_offset
    probe_j = j - eta_offset
    factor = variation_values[:, :, probe_i, probe_j]
    output = np.zeros(
        (factor.shape[0], *reference_values.shape),
        dtype=np.result_type(factor, reference_values),
    )
    diagonal = np.arange(factor.shape[0])
    output[diagonal, diagonal, :, probe_i, probe_j] = (
        factor * reference_values[:, :, probe_i, probe_j] * var_size
    )
    return output


variation_values = {
    "hlt": hlt_eff_var_values,
    "hlt_dtdt": hlt_eff_var_dtdt_values,
    "id": id_eff_var_values,
    "id_dtdt": id_eff_var_dtdt_values,
    "iso": iso_eff_var_values,
    "iso_dtdt": iso_eff_var_dtdt_values,
}
reference_values = {
    "h3": h3_mc_eff.values(),
    "h2": h2_mc_eff.values(),
    "h1": h1_mc_eff.values(),
    "h0": h0_mc_eff.values(),
}
tag_axes = ("time", "mll", "pt_probe", "eta_probe", "pt_tag", "eta_tag")
reference_tag_values = {
    "h3": iso_mc.project(*tag_axes).values(),
    "h2": dtdt_mc.project(*tag_axes).values(),
    "h1": dtst_mc.project(*tag_axes).values(),
    "h0": stst_mc.project(*tag_axes).values(),
}


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

### so at this point i have already selected the mass bin, need to iterate over pt, eta, time

start = time.time()  ############### EFFICIENCY LOOP ##################


def add_batched_systematic(
    values, reference_hist, name, channel, group, as_difference=False, add_poi=False
):
    batched = _batched_hist(values, reference_hist)
    writer.add_systematic(
        batched[{"mll": mass_bin}],
        name,
        "Zmumu",
        channel,
        constrained=False,
        groups=[group],
        as_difference=as_difference,
    )
    if add_poi:
        # values = h3 + variation without the pt_probe flow bins, while h3_poi is
        # projected including them: add only the variation to the h3_poi nominal
        poi = hist.Hist(
            _batched_gen_axis(),
            *h3_poi.axes,
            data=h3_poi.values()[None, ...]
            + (values - reference_values["h3"][None, ...]).sum(
                axis=-2
            ),  # sum over pt_probe
        )
        for low, poi_channel in ((False, "ch_iso_poi_low"), (True, "ch_iso_poi_high")):
            writer.add_systematic(
                remove_bins(poi, "mll", mass_bin + 1, low),
                name,
                "Zmumu",
                poi_channel,
                constrained=False,
                groups=[group],
            )


for i in range(nbins_pt):
    print(f"pt bin: {i}")
    for j in range(nbins_eta):
        print(f"eta_bin: {j}")

        if i > 0:
            writer.add_systematic(
                hist.Hist(
                    *h2_mc_eff.axes,
                    data=_normalization_array(
                        reference_values["h2"], i, j, h2_pt_offset, h2_eta_offset
                    ),
                )[{"mll": mass_bin}],
                f"n_pt{i}_eta{j}",
                "Zmumu",
                "ch_dtdt_eff",
                constrained=False,
                groups=["nz"],
                as_difference=True,
            )

        for reference_key, channel in (
            ("h3", "ch_iso_eff"),
            ("h1", "ch_dtst_eff"),
            ("h0", "ch_stst_eff"),
        ):
            writer.add_systematic(
                hist.Hist(
                    *{"h3": h3_mc_eff, "h1": h1_mc_eff, "h0": h0_mc_eff}[
                        reference_key
                    ].axes,
                    data=_normalization_array(reference_values[reference_key], i, j),
                )[{"mll": mass_bin}],
                f"n_pt{i}_eta{j}",
                "Zmumu",
                channel,
                constrained=False,
                groups=["nz"],
                as_difference=True,
            )

        if i > 0:
            hlt_probe_h3 = _probe_variation_array(
                variation_values["hlt"], reference_values["h3"], i, j
            )
            hlt_tag_h3 = _tag_variation_array(
                variation_values["hlt"], reference_tag_values["h3"], i, j
            )
            add_batched_systematic(
                hlt_tag_h3 + hlt_probe_h3 + reference_values["h3"][None, ...],
                h3_mc_eff,
                f"hlt_pt{i}_eta{j}",
                "ch_iso_eff",
                "eff_trig",
                add_poi=True,
            )

            hlt_tag_h2 = _tag_variation_array(
                variation_values["hlt_dtdt"],
                reference_tag_values["h2"],
                i,
                j,
                h2=True,
                pt_offset=h2_pt_offset,
                eta_offset=h2_eta_offset,
            )
            hlt_probe_h2 = _probe_variation_array(
                variation_values["hlt_dtdt"],
                reference_values["h2"],
                i,
                j,
                pt_offset=h2_pt_offset,
                eta_offset=h2_eta_offset,
            )
            add_batched_systematic(
                2.0 * (hlt_tag_h2 + hlt_probe_h2),
                h2_mc_eff,
                f"hlt_pt{i}_eta{j}",
                "ch_dtdt_eff",
                "eff_trig",
                as_difference=True,
            )
            hlt_tag_h1 = _tag_variation_array(
                variation_values["hlt"], reference_tag_values["h1"], i, j
            )
            hlt_probe_h1 = _probe_variation_array(
                variation_values["hlt"], reference_values["h1"], i, j
            )
            add_batched_systematic(
                2.0 * (hlt_tag_h1 - hlt_probe_h1),
                h1_mc_eff,
                f"hlt_pt{i}_eta{j}",
                "ch_dtst_eff",
                "eff_trig",
                as_difference=True,
            )
            hlt_tag_h0 = _tag_variation_array(
                variation_values["hlt"], reference_tag_values["h0"], i, j
            )
            add_batched_systematic(
                2.0 * hlt_tag_h0,
                h0_mc_eff,
                f"hlt_pt{i}_eta{j}",
                "ch_stst_eff",
                "eff_trig",
                as_difference=True,
            )

            id_tag_h2 = _tag_variation_array(
                variation_values["id_dtdt"],
                reference_tag_values["h2"],
                i,
                j,
                h2=True,
                pt_offset=h2_pt_offset,
                eta_offset=h2_eta_offset,
            )
            id_probe_h2 = _probe_variation_array(
                variation_values["id_dtdt"],
                reference_values["h2"],
                i,
                j,
                pt_offset=h2_pt_offset,
                eta_offset=h2_eta_offset,
            )
            add_batched_systematic(
                id_tag_h2 + id_probe_h2,
                h2_mc_eff,
                f"id_pt{i}_eta{j}",
                "ch_dtdt_eff",
                "eff_id",
                as_difference=True,
            )

            iso_tag_h2 = _tag_variation_array(
                variation_values["iso_dtdt"],
                reference_tag_values["h2"],
                i,
                j,
                h2=True,
                pt_offset=h2_pt_offset,
                eta_offset=h2_eta_offset,
            )
            iso_probe_h2 = _probe_variation_array(
                variation_values["iso_dtdt"],
                reference_values["h2"],
                i,
                j,
                pt_offset=h2_pt_offset,
                eta_offset=h2_eta_offset,
            )
            add_batched_systematic(
                2.0 * (iso_tag_h2 - iso_probe_h2),
                h2_mc_eff,
                f"iso_pt{i}_eta{j}",
                "ch_dtdt_eff",
                "eff_iso",
                as_difference=True,
            )

        iso_tag_h3 = _tag_variation_array(
            variation_values["iso"], reference_tag_values["h3"], i, j
        )
        iso_probe_h3 = _probe_variation_array(
            variation_values["iso"], reference_values["h3"], i, j
        )
        add_batched_systematic(
            iso_tag_h3 + iso_probe_h3 + reference_values["h3"][None, ...],
            h3_mc_eff,
            f"iso_pt{i}_eta{j}",
            "ch_iso_eff",
            "eff_iso",
            add_poi=True,
        )
        iso_tag_h1 = _tag_variation_array(
            variation_values["iso"], reference_tag_values["h1"], i, j
        )
        iso_probe_h1 = _probe_variation_array(
            variation_values["iso"], reference_values["h1"], i, j
        )
        add_batched_systematic(
            2.0 * (iso_tag_h1 - iso_probe_h1),
            h1_mc_eff,
            f"iso_pt{i}_eta{j}",
            "ch_dtst_eff",
            "eff_iso",
            as_difference=True,
        )
        iso_tag_h0 = _tag_variation_array(
            variation_values["iso"], reference_tag_values["h0"], i, j
        )
        add_batched_systematic(
            2.0 * iso_tag_h0,
            h0_mc_eff,
            f"iso_pt{i}_eta{j}",
            "ch_stst_eff",
            "eff_iso",
            as_difference=True,
        )

        id_tag_h3 = _tag_variation_array(
            variation_values["id"], reference_tag_values["h3"], i, j
        )
        id_probe_h3 = _probe_variation_array(
            variation_values["id"], reference_values["h3"], i, j
        )
        add_batched_systematic(
            0.5 * (id_tag_h3 + id_probe_h3) + reference_values["h3"][None, ...],
            h3_mc_eff,
            f"id_pt{i}_eta{j}",
            "ch_iso_eff",
            "eff_id",
            add_poi=True,
        )
        id_tag_h1 = _tag_variation_array(
            variation_values["id"], reference_tag_values["h1"], i, j
        )
        id_probe_h1 = _probe_variation_array(
            variation_values["id"], reference_values["h1"], i, j
        )
        add_batched_systematic(
            id_tag_h1 + id_probe_h1,
            h1_mc_eff,
            f"id_pt{i}_eta{j}",
            "ch_dtst_eff",
            "eff_id",
            as_difference=True,
        )
        id_tag_h0 = _tag_variation_array(
            variation_values["id"], reference_tag_values["h0"], i, j
        )
        id_probe_h0 = _probe_variation_array(
            variation_values["id"], reference_values["h0"], i, j
        )
        add_batched_systematic(
            id_tag_h0 - id_probe_h0,
            h0_mc_eff,
            f"id_pt{i}_eta{j}",
            "ch_stst_eff",
            "eff_id",
            as_difference=True,
        )


for i in range(len(background_syst_names)):
    proc_name = background_proc[i]
    if proc_name == "Zmumu fail gen":
        fgen = True
    else:
        fgen = False
    print("proc_name: %s" % proc_name)
    background_syst(
        writer,
        results,
        background_syst_names[i],
        time_proj_low,
        lumi_scaling,
        [lumi_scaling_h, lumi_scaling_bg],
        proc_name,
        f"bkg_{proc_name}",
        fail_gen=fgen,
    )

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

# HFOC cross detector
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

#### RAMSES cross detector
luminometer_syst(
    writer,
    "ramses",
    iso_ramses,
    dtdt_ramses,
    dtst_ramses,
    stst_ramses,
    "stability",
)

## RAMSES linearity
luminometer_syst(
    writer,
    "ramses",
    iso_sbil_ramses,
    dtdt_sbil_ramses,
    dtst_sbil_ramses,
    stst_sbil_ramses,
    "linearity",
)


delta_t = (time.time() - start) / 3600
print(f"time elapsed: {delta_t} hrs")
outfilename = "liv_rapidity"
if args.randTime:
    outfilename += "_randTime"
if args.sameSignMuon:
    outfilename += "_sameSignMuons"
writer.write(outfolder="./", outfilename=outfilename)
