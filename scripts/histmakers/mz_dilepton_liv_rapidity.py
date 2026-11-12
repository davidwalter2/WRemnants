import os
from datetime import datetime

import hist
import numpy as np

import narf
from narf.lumitools import (
    make_brilcalc_filter_helper,
    make_brilcalc_helper,
    make_lumihelper,
)
from wremnants.production import generator_level_definitions as theory_tools
from wremnants.production import muon_efficiencies_liv, muon_prefiring, muon_selections
from wremnants.production import systematics as syst_tools
from wremnants.production.datasets.dataset_tools import getDatasets
from wremnants.production.histmaker_tools import (
    aggregate_groups,
    scale_to_data,
    write_analysis_output,
)
from wremnants.utilities import common, parsing
from wums import logging

low_pt_cutoff = 15
trigger_pt_cutoff = 25

analysis_label = common.analysis_label(os.path.basename(__file__))
parser, initargs = parsing.common_parser(analysis_label)

parser.add_argument(
    "--randTime",
    type=lambda s: s.lower() in ("true", "1", "yes"),
    default=False,
    help="assign times by event number",
)
parser.add_argument(
    "--sameSignMuon",
    type=lambda s: s.lower() in ("true", "1", "yes"),
    default=False,
    help="whether to select muons with the same sign in events",
)
parser.add_argument(
    "--flipEventNumberSplitting",
    action="store_true",
    help="Flip even with odd event numbers to consider the positive or negative muon as the W-like muon",
)


def make_random_timehelper(filename):
    def to_time(x):
        faketime = int(x[-3:-1]) % 24
        return faketime

    return make_brilcalc_helper(filename, idx=5, action=to_time)


def make_timehelper(filename):
    def to_time(x):
        timestamp = datetime.strptime(x, "%m/%d/%y %H:%M:%S")
        ## calculate julian date since 01.01.2000
        jd = (
            367 * timestamp.year
            - np.floor(7 * (timestamp.year + np.floor((timestamp.month + 9) / 12)) / 4)
            + np.floor(275 * timestamp.month / 9)
            + timestamp.day
            - 730531.5
            + (timestamp.hour + timestamp.minute / 60 + timestamp.second / 3600) / 24
        )
        # calculate greenwich mst, IAU 1982 with T in julian centuries since J2000
        T = jd / 36525
        gmst = (
            67310.54841
            + (876600 * 3600 + 8640184.812866) * T
            + 0.093104 * T**2
            - 6.2e-6 * T**3
        ) % 86400
        gmst /= 3600

        # local sidereal time at CMS: add the east longitude of CMS (6.08 deg = 0.1061 rad)
        lst = (gmst + np.degrees(0.1061) / 15) % 24
        return lst

    return make_brilcalc_helper(filename, idx=2, action=to_time)


def mass_extraction(dataframe, name, root_dataype, filter_name):
    condition = (
        lambda x, idx, f=filter_name: f"Sum({f}) > {idx} ? ROOT::Math::PtEtaPhiMVector({x}_pt[{f}][{idx}], {x}_eta[{f}][{idx}], {x}_phi[{f}][{idx}], wrem::muon_mass) : ROOT::Math::PtEtaPhiMVector(0,0,0,0)"
    )

    new_df = dataframe.Define(f"{name}mu_mom4", condition(f"{root_dataype}", 0))
    new_df = new_df.Define(f"{name}smu_mom4", condition(f"{root_dataype}", 1))

    new_df = new_df.Define(
        f"{name}ll_mom4",
        f"ROOT::Math::PxPyPzEVector({name}mu_mom4)+ROOT::Math::PxPyPzEVector({name}smu_mom4)",
    )
    new_df = new_df.Define(f"{name}mll", f"{name}ll_mom4.mass()")
    new_df = new_df.Define(
        f"{name}pass",
        f"{name}mll >= {low_pt_cutoff} && {name}mll <= 120 && Sum({filter_name})==2",
    )

    return new_df


def luminometer_filter(
    df,
    lumi_name,
    filter_helper,
    helper,
    physics_filter_helper,
):
    ## in the specific luminometer, what are all the luminosity blocks
    df = df.Define(f"blocks_in_{lumi_name}", filter_helper, ["run", "luminosityBlock"])

    df = df.Define(
        "blocks_in_physics", physics_filter_helper, ["run", "luminosityBlock"]
    )
    ### filters to only look at events that were recorded by that specific luminometer. allows any luminometer to be be in that run:fill in PHYSICS
    df_filtered = df.Filter(
        f"(blocks_in_{lumi_name} && blocks_in_physics) && (blocks_in_{lumi_name} == blocks_in_physics)"
    )

    ## instantaneous luminosity read by luminometer
    df_filtered = df_filtered.Define(
        f"lumi_in_{lumi_name}", helper, ["run", "luminosityBlock"]
    )

    df_filtered_hist = df_filtered.HistoBoost(
        f"lumi_{lumi_name}", [axis_date], ["time", f"lumi_in_{lumi_name}"]
    )

    # instantaneous luminosity in physics
    df_filtered_hist_nominal = df_filtered.HistoBoost(
        f"lumi_physics_{lumi_name}", [axis_date], ["time", "lumival"]
    )

    ### take the luminosity measured from this detector and divide it by the fill
    df_filtered = df_filtered.Define(
        f"sbilval_{lumi_name}", f"lumi_in_{lumi_name}/fill_count*1/24"
    )
    ## take the nominal value of the luminosity and divide it by fill
    df_filtered = df_filtered.Define(
        f"sbilval_physics_{lumi_name}", "lumival/fill_count*1/24"
    )

    ### histogram
    df_count_hist = df_filtered.HistoBoost(f"count_{lumi_name}", [axis_date], ["time"])

    df_filtered_hist_sbil = df_filtered.HistoBoost(
        f"sbil_{lumi_name}", [axis_date], ["time", f"sbilval_{lumi_name}"]
    )

    if lumi_name == "pcc":
        return (
            df_filtered_hist,
            df_filtered_hist_nominal,
            df_count_hist,
            df_filtered_hist_sbil,
        )
    else:
        return df_filtered_hist, df_filtered_hist_nominal, df_count_hist


def make_prefire_hists(df, results, name, mass_name="goodMed"):
    df_H = df.Define(
        "weight_newMuonPrefiringSF_H",
        muon_prefiring_helper_H,
        ["Muon_eta", "Muon_pt", "Muon_phi", "Muon_charge", "Muon_tightId"],
    )
    df_BG = df.Define(
        "weight_newMuonPrefiringSF_BG",
        muon_prefiring_helper_BG,
        ["Muon_eta", "Muon_pt", "Muon_phi", "Muon_charge", "Muon_tightId"],
    )

    df_BG = df_BG.Redefine("weight", "weight * weight_newMuonPrefiringSF_BG")
    df_H = df_H.Redefine("weight", "weight * weight_newMuonPrefiringSF_H")

    if mass_name == "goodMed":
        ### technically low and high are the same now,
        a1 = axis_pt_probe
        a2 = axis_pt_tag
        pt_1 = "pt_probe"
        pt_2 = "pt_tag"
        eta_1 = "eta_probe"
        eta_2 = "eta_tag"

        h_weights = df_H.HistoBoost(
            f"{name}_H",
            [
                axis_mll,
                a1,
                axis_eta,
                a2,
                axis_eta_tag,
            ],
            [
                f"{mass_name}_mll",
                pt_1,
                eta_1,
                pt_2,
                eta_2,
                "weight",
            ],
        )
        bg_weights = df_BG.HistoBoost(
            f"{name}_BG",
            [
                axis_mll,
                a1,
                axis_eta,
                a2,
                axis_eta_tag,
            ],
            [
                f"{mass_name}_mll",
                pt_1,
                eta_1,
                pt_2,
                eta_2,
                "weight",
            ],
        )
        syst_tools.add_L1Prefire_unc_hists(
            results,
            df_BG,
            [
                axis_mll,
                a1,
                axis_eta,
                a2,
                axis_eta_tag,
            ],
            [
                f"{mass_name}_mll",
                pt_1,
                eta_1,
                pt_2,
                eta_2,
            ],
            helper_stat=muon_prefiring_helper_stat_BG,
            helper_syst=muon_prefiring_helper_syst_BG,
            storage_type=hist.storage.Double(),
            base_name=name + "_BG",
            weight="weight",
        )

        syst_tools.add_L1Prefire_unc_hists(
            results,
            df_H,
            [
                axis_mll,
                a1,
                axis_eta,
                a2,
                axis_eta_tag,
            ],
            [
                f"{mass_name}_mll",
                pt_1,
                eta_1,
                pt_2,
                eta_2,
            ],
            helper_stat=muon_prefiring_helper_stat_H,
            helper_syst=muon_prefiring_helper_syst_H,
            storage_type=hist.storage.Double(),
            base_name=name + "_H",
            weight="weight",
        )
    else:
        a1 = axis_pt_probe
        pt_1 = "Muon_pt"
        eta_1 = "Muon_eta"

        h_weights = df_H.HistoBoost(
            f"{name}_H",
            [
                axis_mll,
                a1,
                axis_eta,
            ],
            [
                f"{mass_name}_mll",
                pt_1,
                eta_1,
                "weight",
            ],
        )
        bg_weights = df_BG.HistoBoost(
            f"{name}_BG",
            [
                axis_mll,
                a1,
                axis_eta,
            ],
            [
                f"{mass_name}_mll",
                pt_1,
                eta_1,
                "weight",
            ],
        )

        syst_tools.add_L1Prefire_unc_hists(
            results,
            df_BG,
            [
                axis_mll,
                a1,
                axis_eta,
            ],
            [
                f"{mass_name}_mll",
                pt_1,
                eta_1,
            ],
            helper_stat=muon_prefiring_helper_stat_BG,
            helper_syst=muon_prefiring_helper_syst_BG,
            storage_type=hist.storage.Double(),
            base_name=name + "_BG",
            weight="weight",
        )

        syst_tools.add_L1Prefire_unc_hists(
            results,
            df_H,
            [
                axis_mll,
                a1,
                axis_eta,
            ],
            [
                f"{mass_name}_mll",
                pt_1,
                eta_1,
            ],
            helper_stat=muon_prefiring_helper_stat_H,
            helper_syst=muon_prefiring_helper_syst_H,
            storage_type=hist.storage.Double(),
            base_name=name + "_H",
            weight="weight",
        )
    results.append(h_weights)
    results.append(bg_weights)
    return bg_weights, h_weights


def true_efficiencies(df):
    df = df.Define(
        "reco_scalefactor_weight",
        muon_reco_efficiency_helper,
        ["Muon_eta", "Muon_pt"],
    )
    df = df.Define(
        "tracking_scalefactor_weight",
        muon_tracking_efficiency_helper,
        ["Muon_standaloneEta", "Muon_standalonePt"],
    )
    df = df.Redefine(
        "weight", "weight * reco_scalefactor_weight* tracking_scalefactor_weight"
    )

    df = df.Define("PosMuon_Global", "Muon_isPositive && Muon_isGoodGlobal")
    df = df.Define("PosMuon_ID", "PosMuon_Global && (Muon_isGoodMedium == 1)")
    df = df.Define("PosMuon_trig", "PosMuon_ID && Muon_isGoodTrigger")

    df = df.Define(
        "PosMuon_iso",
        f"Muon_passIso == 1 && ((PosMuon_trig && Muon_pt >= {trigger_pt_cutoff}) || (Muon_pt < {trigger_pt_cutoff} && PosMuon_ID))",
    )

    df = df.Define("Global_pt", "Muon_pt[PosMuon_Global]")
    df = df.Define("Global_eta", "Muon_eta[PosMuon_Global]")

    df = df.Define("ID_pt", "Muon_pt[PosMuon_ID]")
    df = df.Define("ID_eta", "Muon_eta[PosMuon_ID]")

    df = df.Define("Trig_pt", "Muon_pt[PosMuon_trig]")
    df = df.Define("Trig_eta", "Muon_eta[PosMuon_trig]")

    df = df.Define("Iso_pt", "Muon_pt[PosMuon_iso]")
    df = df.Define("Iso_eta", "Muon_eta[PosMuon_iso]")
    return df


axis_date = hist.axis.Regular(24, 0, 24, name="time", overflow=False, underflow=False)

axis_sbil = hist.axis.Regular(
    24, 9e-7, 3e-8, name="sbil", overflow=False, underflow=False
)


axis_eta = hist.axis.Variable([-2.4, -1.4, -0.7, 0, 0.7, 1.4, 2.4], name="eta_probe")
axis_eta_tag = hist.axis.Variable([-2.4, -1.4, -0.7, 0, 0.7, 1.4, 2.4], name="eta_tag")


axis_pt_tag = hist.axis.Variable(
    [
        low_pt_cutoff,
        trigger_pt_cutoff,
        28,
        30,
        32,
        34,
        36,
        38,
        40,
        42,
        47,
        55,
        60,
        65,
        80,
    ],
    name="pt_tag",
)

axis_pt_tag_copy = hist.axis.Variable(
    [
        low_pt_cutoff,
        trigger_pt_cutoff,
        28,
        30,
        32,
        34,
        36,
        38,
        40,
        42,
        47,
        55,
        60,
        65,
        80,
    ],
    name="pt_tag",
)

axis_pt_probe = hist.axis.Variable(
    [
        low_pt_cutoff,
        trigger_pt_cutoff,
        28,
        30,
        32,
        34,
        36,
        38,
        40,
        42,
        47,
        55,
        60,
        65,
        80,
    ],
    name="pt_probe",
)

axis_pt_probe_copy = hist.axis.Variable(
    [
        low_pt_cutoff,
        trigger_pt_cutoff,
        28,
        30,
        32,
        34,
        36,
        38,
        40,
        42,
        47,
        55,
        60,
        65,
        80,
    ],
    name="pt_probe",
)

## REMEMBER TO SWITCH BACK TO THIS
axis_mll = hist.axis.Variable(
    [15, 30, 40, 45, 50, 55, 60, 65, 70, 76, 106, 110, 115, 120], name="mll"
)
axis_mll_copy = hist.axis.Variable(
    [15, 30, 40, 45, 50, 55, 60, 65, 70, 76, 106, 110, 115, 120], name="goodMed_mll"
)

axis_mVgen = hist.axis.Variable(
    [15, 30, 40, 45, 50, 55, 60, 65, 70, 76, 106, 110, 115, 120],
    name="genMass",
    flow=False,
)
axis_mreco = hist.axis.Variable(
    [15, 30, 40, 45, 50, 55, 60, 65, 70, 76, 106, 110, 115, 120],
    name="recoMass",
    flow=False,
)

axis_weight = hist.axis.Regular(50, 0.5, 1, name="weight")

args = parser.parse_args()
logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)
era = args.era
calib_filepaths = common.calib_filepaths
lumi_files_path = f"{common.base_dir}/wremnants/production/datasets"  #### THIS IS A REALLY DUMB WAY TO DO THIS
# hoping this can go up top
lumicsv = f"{lumi_files_path}/bylsoutput_nBunches.csv"
hfoc_csv = f"{lumi_files_path}/bylsoutput_nBunches_HFOC.csv"
pcc_csv = f"{lumi_files_path}/bylsoutput_nBunches_PCC.csv"
ramses_csv = f"{lumi_files_path}/bylsoutput_nBunches_RAMSES.csv"

brilcalc_helper = make_timehelper(lumicsv)
lumi_no_time = make_lumihelper(lumicsv)  # post_vfp
lumi_bunch_helper = make_brilcalc_helper(lumicsv, idx=9, action=float)

randomTime_helper = make_random_timehelper(lumicsv)

### lumihelpers return the instantaneous luminosity for a given run:fill
hfoc_helper = make_lumihelper(hfoc_csv)
pcc_helper = make_lumihelper(pcc_csv)
ramses_helper = make_lumihelper(ramses_csv)

# filter helpers return run:fill for a certain csv so that other csvs can be filtered by that
hfoc_filter_helper = make_brilcalc_filter_helper(hfoc_csv)
pcc_filter_helper = make_brilcalc_filter_helper(pcc_csv)
ramses_filter_helper = make_brilcalc_filter_helper(ramses_csv)
physics_filter_helper = make_brilcalc_filter_helper(lumicsv)


## make the two sets of prefiring helpers for each port of the data
(
    muon_prefiring_helper_BG,
    muon_prefiring_helper_stat_BG,
    muon_prefiring_helper_syst_BG,
) = muon_prefiring.make_muon_prefiring_helpers(era="2016BG")

muon_prefiring_helper_H, muon_prefiring_helper_stat_H, muon_prefiring_helper_syst_H = (
    muon_prefiring.make_muon_prefiring_helpers(era="2016H")
)


muon_reco_efficiency_helper = muon_efficiencies_liv.make_muon_reco_efficiency_helper(
    era="2016BG",
)
muon_tracking_efficiency_helper = (
    muon_efficiencies_liv.make_muon_tracking_efficiency_helper(
        era="2016BG",
    )
)

datasets = getDatasets(
    maxFiles=args.maxFiles,
    filt=args.filterProcs,
    excl=args.excludeProcs,
    nanoVersion="v9",
    base_path=args.dataPath,
    extended="msht20an3lo" not in args.pdfs,
    era=era,
)


########################################################
def build_graph_lumi(df, dataset):
    ## get the physics values

    df = df.Define("time", brilcalc_helper, ["run", "luminosityBlock"])
    hist_lumi_nom = df.HistoBoost("lumi_nom", [axis_date], ["time", "lumival"])
    df = df.Define("fill_count", lumi_bunch_helper, ["run", "luminosityBlock"])

    df_H = df.Filter("run >= 281613")
    hist_lumi_post = df_H.HistoBoost("lumi_post", [axis_date], ["time", "lumival"])
    df_B = df.Filter("run < 281613")
    hist_lumi_pre = df_B.HistoBoost("lumi_pre", [axis_date], ["time", "lumival"])

    hist_hfoc_filtered, hist_nominal_in_hfoc_filtered, hist_hfoc_count = (
        luminometer_filter(
            df,
            "hfoc",
            hfoc_filter_helper,
            hfoc_helper,
            physics_filter_helper,
        )
    )
    hist_pcc_filtered, hist_nominal_in_pcc_filtered, hist_pcc_count, hist_pcc_sbil = (
        luminometer_filter(
            df,
            "pcc",
            pcc_filter_helper,
            pcc_helper,
            physics_filter_helper,
        )
    )
    hist_ramses_filtered, hist_nominal_in_ramses_filtered, hist_ramses_count = (
        luminometer_filter(
            df,
            "ramses",
            ramses_filter_helper,
            ramses_helper,
            physics_filter_helper,
        )
    )

    results = [
        hist_lumi_nom,
        hist_hfoc_filtered,
        hist_pcc_filtered,
        hist_ramses_filtered,
        hist_nominal_in_hfoc_filtered,
        hist_nominal_in_pcc_filtered,
        hist_nominal_in_ramses_filtered,
        hist_pcc_sbil,
        hist_hfoc_count,
        hist_pcc_count,
        hist_ramses_count,
        hist_lumi_pre,
        hist_lumi_post,
    ]
    return results


def build_graph(df, dataset):
    isoBranch = muon_selections.getIsoBranch("iso04vtxAgn")
    logger.info(f"fomrbuild graph for dataset: {dataset.name}")
    results = []

    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
        if not args.randTime:
            df = df.Define("time", brilcalc_helper, ["run", "luminosityBlock"])
        else:
            df = df.Define("time", "rand()%24")

        hist_time = df.HistoBoost("time", [axis_date], ["time"])
    else:
        ### oh do I need this?
        # weight_expr += "*weight_fullMuonSF_withTrackingReco"
        df = df.Define("weight", "std::copysign(1.0, genWeight)")
    weightsum = df.SumAndCount("weight")

    ###### NEED TO GENERATE THE RIGHT VARIABLES ####    )
    df = df.Define(
        "isEvenEvent", f"event % 2 {'!=' if args.flipEventNumberSplitting else '=='} 0"
    )

    ### same as mW, allow for lower pt
    df = df.Define(
        "Muon_isGoodGlobal",
        f" Muon_isGlobal && Muon_highPurity && Muon_standaloneNumberOfValidHits > 0 && Muon_standalonePt > {low_pt_cutoff} &&  wrem::vectDeltaR2(Muon_standaloneEta, Muon_standalonePhi, Muon_eta, Muon_phi) < 0.09 && Muon_pt >= {low_pt_cutoff} && abs(Muon_eta) <= 2.4 && Muon_charge != -99 && abs(Muon_dxybs) < 0.05 && Muon_looseId",
    )

    df = df.Define(
        "Muon_isGoodMedium",
        f"Muon_isGoodGlobal && Muon_mediumId ",
    )

    df = df.Define(
        "goodTrigObjs",
        f"wrem::goodMuonTriggerCandidate<wrem::Era::Era_2016PostVFP>(TrigObj_id,TrigObj_filterBits)",
    )

    df = df.Define(
        "Muon_isGoodTrigger",
        f"Muon_pt>={trigger_pt_cutoff} && Muon_isGoodMedium && wrem::hasTriggerMatch(Muon_eta,Muon_phi,TrigObj_eta[goodTrigObjs],TrigObj_phi[goodTrigObjs])",
    )

    df = df.Define("Muon_passIso", f"({isoBranch} < {args.isolationThreshold})")

    df = df.Define(
        "Muon_passIsoTrig", "(Muon_isGoodTrigger == 1) && (Muon_passIso == 1)"
    )

    #### all the cut types are fully defined above. now actually doing the cuts

    ### I have this set aside for the true efficiency comparison
    df_positive = df.Filter("Sum(Muon_isGoodGlobal) > 0 ")
    df_positive = df_positive.Define(
        "Muon_isPositive", "Muon_charge == 1 || Muon_charge == -1 "
    )

    #### always want events with only two muons
    df = df.Filter("Sum(Muon_isGoodGlobal) == 2")

    ### require that they are opposite charges
    if not args.sameSignMuon:
        df = df.Filter(
            "Muon_charge[Muon_isGoodGlobal][0] != Muon_charge[Muon_isGoodGlobal][1]"
        )
    else:
        df = df.Filter(
            "Muon_charge[Muon_isGoodGlobal][0] == Muon_charge[Muon_isGoodGlobal][1]"
        )
    #### filter to ensure that at least one muon passes isolation
    df = df.Filter(
        "(Muon_passIsoTrig[Muon_isGoodGlobal][0] == 1) || (Muon_passIsoTrig[Muon_isGoodGlobal][1] == 1)"
    )
    ### define which one is the probe. if both pass all cuts, decide according to which one is even
    # if muon0 passes trigger & tight id (if muon1 passes trigger and tight id: then randomly select 0 vs 1, if muon1 fails, make it the probe), if muon0 fails make it 0
    df = df.Define(
        "mu_probe",
        "Muon_passIsoTrig[Muon_isGoodGlobal][0] == 1 ? (Muon_passIsoTrig[Muon_isGoodGlobal][1] == 1 ? isEvenEvent: 1) : 0",
    )

    ### other muon should be the tag
    df = df.Define("mu_tag", "mu_probe == 0 ? 1: 0")

    ### getting masses and then assigning them appropriately
    df = mass_extraction(df, "goodMed_", "Muon", "Muon_isGoodGlobal")

    df = df.Define(
        "pt_probe", "mu_probe == 0 ? goodMed_mu_mom4.pt() : goodMed_smu_mom4.pt()"
    )
    df = df.Define(
        "eta_probe",
        "mu_probe == 0 ? goodMed_mu_mom4.eta() : goodMed_smu_mom4.eta()",
    )

    df = df.Define(
        "pt_tag", "mu_probe == 0 ? goodMed_smu_mom4.pt() : goodMed_mu_mom4.pt()"
    )
    df = df.Define(
        "eta_tag", "mu_probe == 0 ? goodMed_smu_mom4.eta():goodMed_mu_mom4.eta()"
    )

    if not dataset.is_data:

        df_positive = true_efficiencies(df_positive)
        df = df.Define(
            "reco_scalefactor_weight",
            muon_reco_efficiency_helper,
            ["Muon_eta", "Muon_pt"],
        )
        df = df.Define(
            "tracking_scalefactor_weight",
            muon_tracking_efficiency_helper,
            ["Muon_standaloneEta", "Muon_standalonePt"],
        )

        ### need to correct for these
        df = df.Redefine(
            "weight", "weight * reco_scalefactor_weight * tracking_scalefactor_weight"
        )
        df = theory_tools.define_postfsr_vars(df)
        df = theory_tools.define_prefsr_vars(df)

        df_positive = theory_tools.define_postfsr_vars(df_positive)
        ### positive stuff is just for comparison with mw efficincies

        response_matrix = df.HistoBoost(
            "mass_response",
            [
                axis_mVgen,
                axis_mreco,
            ],
            [
                "massVgen",
                "goodMed_mll",
                "weight",
            ],
        )

        pos_global = df_positive.HistoBoost(
            "pos_global",
            [
                axis_pt_probe,
                axis_eta,
            ],
            [
                "Global_pt",
                "Global_eta",
                "weight",
            ],
        )

        pos_ID = df_positive.HistoBoost(
            "pos_ID",
            [
                axis_pt_probe,
                axis_eta,
            ],
            [
                "ID_pt",
                "ID_eta",
                "weight",
            ],
        )

        pos_trig = df_positive.HistoBoost(
            "pos_trig",
            [
                axis_pt_probe,
                axis_eta,
            ],
            [
                "Trig_pt",
                "Trig_eta",
                "weight",
            ],
        )

        pos_iso = df_positive.HistoBoost(
            "pos_iso",
            [
                axis_pt_probe,
                axis_eta,
            ],
            [
                "Iso_pt",
                "Iso_eta",
                "weight",
            ],
        )

        ### back to doing things for the in situ measuremetns
        df = df.Define(
            "postfsrMuons_loose",
            f"postfsrMuons && abs(GenPart_eta) < 2.4 && GenPart_pt > {low_pt_cutoff}",
        )

        df = mass_extraction(df, "gen_", "GenPart", "postfsrMuons_loose")

        # df_loose = df.Filter("gen_pass")
        df_loose = df
        hist_pass_gen = df.HistoBoost(
            "pass_gen",
            [
                axis_mll,
                axis_pt_probe,
                axis_eta,
                axis_pt_tag,
                axis_eta_tag,
            ],
            [
                "goodMed_mll",
                "pt_probe",
                "eta_probe",
                "pt_tag",
                "eta_tag",
                "weight",
            ],
        )

        df_tight = df_loose.Filter(
            "Muon_isGoodMedium[Muon_isGoodGlobal][mu_probe] == 1"
        )
        df_trig = df_tight.Filter(
            "Muon_isGoodTrigger[Muon_isGoodGlobal][mu_probe] == 1"
        )

        ## either the probe muon passes the trigger or it is going to pass iso. or is okay because trigger is above 25 and this is mutually exclusive
        df_iso = df_tight.Filter(
            f"Muon_passIso[Muon_isGoodGlobal][mu_probe] == 1 && ((Muon_isGoodTrigger[Muon_isGoodGlobal][mu_probe] == 1 && Muon_pt[Muon_isGoodGlobal][mu_probe] > {trigger_pt_cutoff}) || Muon_pt[Muon_isGoodGlobal][mu_probe] < {trigger_pt_cutoff})"
        )

        # ### fail generator
        df_loose_fg = df.Filter("!gen_pass")
        df_tight_fg = df_loose_fg.Filter(
            "Muon_isGoodMedium[Muon_isGoodGlobal][mu_probe] == 1"
        )
        df_trig_fg = df_tight_fg.Filter(
            "Muon_isGoodTrigger[Muon_isGoodGlobal][mu_probe] == 1"
        )

        results.append(hist_pass_gen)
        results.append(pos_global)
        results.append(pos_ID)
        results.append(pos_trig)
        results.append(pos_iso)
        results.append(response_matrix)

        make_prefire_hists(df_iso, results, "pass_iso")
        make_prefire_hists(df_trig, results, "dtdt_prpg")
        make_prefire_hists(df_tight, results, "dtst_prpg")
        make_prefire_hists(df_loose, results, "stst_prpg")
        make_prefire_hists(df_trig_fg, results, "dtdt_prfg")
        make_prefire_hists(df_tight_fg, results, "dtst_prfg")
        make_prefire_hists(df_loose_fg, results, "stst_prfg")

    else:  ### this is for real data

        stst = df
        dtst = stst.Filter("Muon_isGoodMedium[Muon_isGoodGlobal][mu_probe] == 1")
        dtdt = dtst.Filter("Muon_isGoodTrigger[Muon_isGoodGlobal][mu_probe] == 1")

        iso = dtst.Filter(
            f"Muon_passIso[Muon_isGoodGlobal][mu_probe] == 1 && ((Muon_isGoodTrigger[Muon_isGoodGlobal][mu_probe] == 1 && Muon_pt[Muon_isGoodGlobal][mu_probe] > {trigger_pt_cutoff}) || Muon_pt[Muon_isGoodGlobal][mu_probe] < {trigger_pt_cutoff})"
        )

        hist_time_proj = df.HistoBoost(
            "time_proj",
            [
                axis_date,
                axis_mll,
                axis_pt_probe,
                axis_eta,
                axis_pt_tag,
                axis_eta_tag,
            ],
            [
                "time",
                "goodMed_mll",
                "pt_probe",
                "eta_probe",
                "pt_tag",
                "eta_tag",
                "weight",
            ],
        )

        hist_time_iso = iso.HistoBoost(
            "time_iso",
            [
                axis_date,
                axis_mll,
                axis_pt_probe,
                axis_eta,
                axis_pt_tag,
                axis_eta_tag,
            ],
            [
                "time",
                "goodMed_mll",
                "pt_probe",
                "eta_probe",
                "pt_tag",
                "eta_tag",
                "weight",
            ],
        )

        hist_time_mll = dtdt.HistoBoost(
            "time_mll",
            [
                axis_date,
                axis_mll,
                axis_pt_probe,
                axis_eta,
                axis_pt_tag,
                axis_eta_tag,
            ],
            [
                "time",
                "goodMed_mll",
                "pt_probe",
                "eta_probe",
                "pt_tag",
                "eta_tag",
                "weight",
            ],
        )
        hist_time_dtst = dtst.HistoBoost(
            "time_dtst",
            [
                axis_date,
                axis_mll,
                axis_pt_probe,
                axis_eta,
                axis_pt_tag,
                axis_eta_tag,
            ],
            [
                "time",
                "goodMed_mll",
                "pt_probe",
                "eta_probe",
                "pt_tag",
                "eta_tag",
                "weight",
            ],
        )
        hist_time_stst = stst.HistoBoost(
            "time_stst",
            [
                axis_date,
                axis_mll,
                axis_pt_probe,
                axis_eta,
                axis_pt_tag,
                axis_eta_tag,
            ],
            [
                "time",
                "goodMed_mll",
                "pt_probe",
                "eta_probe",
                "pt_tag",
                "eta_tag",
                "weight",
            ],
        )

        results.append(hist_time_iso)
        results.append(hist_time_proj)
        results.append(hist_time)
        results.append(hist_time_mll)
        results.append(hist_time_dtst)
        results.append(hist_time_stst)

    return results, weightsum


logger.debug(f"Datasets are {[d.name for d in datasets]}")
resultdict = narf.build_and_run(datasets[::-1], build_graph, build_graph_lumi)

if not args.noScaleToData:
    scale_to_data(resultdict)
    aggregate_groups(datasets, resultdict, args.aggregateGroups)

write_analysis_output(
    resultdict, f"{os.path.basename(__file__).replace('py', 'hdf5')}", args
)
