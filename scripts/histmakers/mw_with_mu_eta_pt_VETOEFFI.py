import os

import numpy as np

from utilities import common, parsing
from wremnants.datasets.datagroups import Datagroups
from wums import logging

analysis_label = Datagroups.analysisLabel(
    os.path.basename(__file__).replace("_VETOEFFI", "")
)
parser, initargs = parsing.common_parser(analysis_label)

import math
import os

import hist
import ROOT

import narf
from wremnants import (
    muon_calibration,
    muon_selections,
    pileup,
    theory_corrections,
    vertex,
)
from wremnants.datasets.dataset_tools import getDatasets
from wremnants.histmaker_tools import (
    aggregate_groups,
    scale_to_data,
    write_analysis_output,
)

parser.add_argument(
    "--vetoGenPartPt",
    type=float,
    default=0.0,
    help="Minimum pT for the postFSR gen muon",
)
parser.add_argument(
    "--makeCorr",
    action="store_true",
    help="Make muon drop correction helper",
)
parser.add_argument(
    "--selectReco",
    action="store_true",
    help="Apply reco muon selection",
)
parser.add_argument(
    "--tnpEfficiency",
    action="store_true",
    help="Compute T&P like efficiency",
)
args = parser.parse_args()
logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

args = parser.parse_args()

thisAnalysis = ROOT.wrem.AnalysisType.Wmass

era = args.era
datasets = getDatasets(
    maxFiles=args.maxFiles,
    filt=args.filterProcs,
    excl=args.excludeProcs,
    nanoVersion="v9",
    base_path=args.dataPath,
    oneMCfileEveryN=args.oneMCfileEveryN,
    extended="msht20an3lo" not in args.pdfs,
    era=era,
)

# # custom template binning
# template_neta = int(args.eta[0])
# template_mineta = args.eta[1]
# template_maxeta = args.eta[2]
# logger.info(
#     f"Eta binning: {template_neta} bins from {template_mineta} to {template_maxeta}"
# )
# template_npt = int(args.pt[0])
# template_minpt = args.pt[1]
# template_maxpt = args.pt[2]
# logger.info(
#     f"Pt binning: {template_npt} bins from {template_minpt} to {template_maxpt}"
# )

veto_minpt = 15
veto_maxpt = 60
veto_npt = int(veto_maxpt - veto_minpt)


veto_maxeta = 2.4
veto_neta = int(10 * veto_maxeta * 2)

# standard regular axes
axis_eta0 = hist.axis.Regular(
    veto_neta,
    -veto_maxeta,
    veto_maxeta,
    name="eta",
    overflow=True,
    underflow=True,
)
axis_pt0 = hist.axis.Regular(
    veto_npt,
    veto_minpt,
    veto_maxpt,
    name="pt",
    overflow=True,
    underflow=True,
)
axis_eta1 = hist.axis.Regular(
    veto_neta,
    -veto_maxeta,
    veto_maxeta,
    name="eta1",
    overflow=True,
    underflow=True,
)
axis_pt1 = hist.axis.Regular(
    veto_npt,
    veto_minpt,
    veto_maxpt,
    name="pt1",
    overflow=True,
    underflow=True,
)

# for the observed contributions
axis_eta = hist.axis.Regular(
    48,
    -2.4,
    2.4,
    name="eta",
    overflow=False,
    underflow=False,
)
axis_pt = hist.axis.Regular(
    30,
    26,
    56,
    name="pt",
    overflow=False,
    underflow=False,
)

axis_charge = common.axis_charge
axis_passVeto = hist.axis.Boolean(name="passVeto")
axis_phi = hist.axis.Regular(25, 0, math.pi, name="phi", overflow=True)

axis_met = hist.axis.Regular(25, 0.0, 100.0, name="met", underflow=False, overflow=True)
axis_mt = hist.axis.Variable(
    (*np.arange(0, 95, 5), 100, 120),
    name="mt",
    underflow=False,
    overflow=True,
)

axis_relIso = hist.axis.Regular(
    100, 0, 1, name="relIso", underflow=False, overflow=True
)
axis_dxy = hist.axis.Regular(100, 0, 0.01, name="dxy", underflow=False, overflow=True)

nominal_axes = [axis_eta, axis_pt, axis_charge, axis_passVeto]
nominal_cols = [
    "postFSRmuon_eta0",
    "postFSRmuon_pt0",
    "postFSRmuon_charge0",
    "passVeto",
]

# sum those groups up in post processing
groups_to_aggregate = args.aggregateGroups

# qcdScaleByHelicity_helper = theory_corrections.make_qcd_uncertainty_helper_by_helicity()

logger.info("Running with no scale factors")

pileup_helper = pileup.make_pileup_helper(era=era)
vertex_helper = vertex.make_vertex_helper(era=era)

calib_filepaths = common.calib_filepaths
closure_filepaths = common.closure_filepaths

(
    mc_jpsi_crctn_helper,
    data_jpsi_crctn_helper,
    jpsi_crctn_MC_unc_helper,
    jpsi_crctn_data_unc_helper,
) = muon_calibration.make_jpsi_crctn_helpers(
    args, calib_filepaths, make_uncertainty_helper=True
)

mc_calibration_helper, data_calibration_helper, calibration_uncertainty_helper = (
    muon_calibration.make_muon_calibration_helpers(args)
)

smearing_helper, smearing_uncertainty_helper = (
    (None, None) if args.noSmearing else muon_calibration.make_muon_smearing_helpers()
)

bias_helper = (
    muon_calibration.make_muon_bias_helpers(args) if args.biasCalibration else None
)

theory_corrs = [*args.theoryCorr, *args.ewTheoryCorr]
procsWithTheoryCorr = [d.name for d in datasets if d.name in common.vprocs]
if len(procsWithTheoryCorr):
    corr_helpers = theory_corrections.load_corr_helpers(
        procsWithTheoryCorr, theory_corrs
    )
else:
    corr_helpers = {}

smearing_weights_procs = []


import h5py
import numpy as np

from utilities.io_tools import input_tools
from wremnants.correctionsTensor_helper import makeCorrectionsTensor
from wums import boostHistHelpers as hh

if not args.makeCorr:

    # filename = "/home/submit/david_w/ceph/Unfolding/results_histmaker/2508011_muonDrop/mw_with_mu_eta_pt_VETOEFFI_scetlib_dyturboCorr.hdf5"
    filename = "/home/submit/david_w/ceph/Unfolding/results_histmaker/2508013_muonDrop/mw_with_mu_eta_pt_VETOEFFI_scetlib_dyturboCorr.hdf5"
    h5file = h5py.File(filename, "r")
    results = input_tools.load_results_h5py(h5file)

    def loadHist(sample, name):
        h = results[sample]["output"][name].get()
        scale = results[sample]["dataset"]["xsec"] / results[sample]["weight_sum"]
        h = hh.scaleHist(h, scale)
        return h

    # low mass
    h1 = loadHist("DYJetsToMuMuMass10to50PostVFP", "drop")

    hp = h1[{"pt": slice(26j, 56j), "eta": slice(-2.4j, 2.4j)}]
    hp_all = hp[
        {"eta1": slice(None, None, hist.sum), "pt1": slice(None, None, hist.sum)}
    ]
    hp_in = hp[
        {
            "eta1": slice(0, hist.overflow, hist.sum),
            "pt1": slice(0, None, hist.sum),
        }
    ]
    # observerd positive muon inside, other ourside
    hp_out = hh.addHists(hp_all, hp_in, scale2=-1)

    hm = h1[{"pt1": slice(26j, 56j), "eta1": slice(-2.4j, 2.4j)}]
    hm_all = hm[{"eta": slice(None, None, hist.sum), "pt": slice(None, None, hist.sum)}]
    hm_in = hm[
        {
            "eta": slice(0, hist.overflow, hist.sum),
            "pt": slice(0, None, hist.sum),
        }
    ]
    # observerd negative muon inside, other ourside
    hm_out = hh.addHists(hm_all, hm_in, scale2=-1)

    # high mass
    h2 = loadHist("ZmumuPostVFP", "drop")

    h2p = h2[{"pt": slice(26j, 56j), "eta": slice(-2.4j, 2.4j)}]
    h2m = h2[{"pt1": slice(26j, 56j), "eta1": slice(-2.4j, 2.4j)}]
    h2p_in = h2p[
        {
            "eta1": slice(0, hist.overflow, hist.sum),
            "pt1": slice(0, None, hist.sum),
        }
    ]
    h2m_in = h2m[
        {
            "eta": slice(0, hist.overflow, hist.sum),
            "pt": slice(0, None, hist.sum),
        }
    ]

    # transfer evetns from high mass to low mass
    hp_prob = hh.divideHists(hp_out.project("eta", "pt"), h2p_in.project("eta", "pt"))
    hm_prob = hh.divideHists(
        hm_out.project("eta1", "pt1"), h2m_in.project("eta1", "pt1")
    )

    axis_dummy = hist.axis.Regular(1, 0, 2, underflow=False, overflow=False, name="var")
    h_drop = hist.Hist(
        *hp_prob.axes, axis_charge, axis_dummy, storage=hist.storage.Double()
    )
    h_drop.values(flow=True)[...] = np.stack(
        [hm_prob.values(flow=True), hp_prob.values(flow=True)], axis=-1
    )[..., np.newaxis]
    drop_helper = makeCorrectionsTensor(h_drop)

    if args.selectReco:
        # for veto efficiency
        eff_helper = {}

        # gen efficiencies
        for key, res in results.items():
            if key not in [d.name for d in datasets]:
                continue

            hp = res["output"]["vetoMuonPass0"].get()
            hm = res["output"]["vetoMuonPass1"].get()

            # transfer factor tensor with two bins: tf = (eff, 1-eff)
            # split by charge
            tfp = (
                hp.values(flow=True)
                / hp[{"passVeto": hist.sum}].values(flow=True)[..., np.newaxis]
            )
            tfm = (
                hm.values(flow=True)
                / hm[{"passVeto": hist.sum}].values(flow=True)[..., np.newaxis]
            )
            h_tf = hist.Hist(axis_charge, *hp.axes, storage=hist.storage.Double())
            h_tf.values(flow=True)[...] = np.stack((tfm, tfp), axis=0)

            ## charge integrated
            # h1 = hh.addHists(hp, hm)

            # tf = (
            #     h1.values(flow=True) / h1[{"passVeto": hist.sum}].values(flow=True)[..., np.newaxis]
            # )

            # h_tf = hist.Hist(*hp.axes, storage=hist.storage.Double())
            # h_tf.values(flow=True)[...] = tf

            eff_helper[key] = makeCorrectionsTensor(h_tf)

        # use efficiencies for Z
        eff_helper["DYlowMass"] = eff_helper["ZmumuPostVFP"]
        eff_helper["QCDmuEnrichPt15PostVFP"] = eff_helper["ZmumuPostVFP"]

        # reco efficiencies
        if args.tnpEfficiency:
            base_dir = f"{common.data_dir}/muonSF/efficiencies_GtoH/"

            def make_efficiency_hist(eff):
                infile_p = f"{base_dir}/mu_{eff}_plus.txt"
                infile_m = f"{base_dir}/mu_{eff}_minus.txt"

                import pandas as pd

                def read_file(infile):
                    with open(infile) as f:
                        var1_name = f.readline().split(":")[1].strip()
                        var2_name = f.readline().split(":")[1].strip()

                    df = pd.read_csv(infile, sep="\t", skiprows=2, engine="python")
                    df.columns = df.columns.str.strip()
                    df = df.apply(pd.to_numeric, errors="raise")

                    var1_edges = np.array(
                        sorted(set(df["var1min"]).union(df["var1max"]))
                    )
                    var2_edges = np.array(
                        sorted(set(df["var2min"]).union(df["var2max"]))
                    )

                    def make_hist(colname):
                        h = hist.Hist(
                            hist.axis.Variable(var1_edges, name=var1_name),
                            hist.axis.Variable(var2_edges, name=var2_name),
                            storage=hist.storage.Double(),
                        )
                        # Fill manually from bins
                        for _, row in df.iterrows():
                            h.fill(
                                **{
                                    var1_name: 0.5 * (row["var1min"] + row["var1max"]),
                                    var2_name: 0.5 * (row["var2min"] + row["var2max"]),
                                },
                                weight=row[colname],
                            )
                        return h

                    h_eff_mc = make_hist("eff mc")
                    # h_eff_data = make_hist("eff data")
                    # h_err_data = make_hist("err data")

                    return h_eff_mc

                hp = read_file(infile_p)
                hm = read_file(infile_m)

                h_eff = hist.Hist(
                    axis_charge,
                    *hp.axes,
                    axis_dummy,
                    storage=hist.storage.Double(),
                    data=np.stack([hp.values(flow=True), hm.values(flow=True)], axis=0)[
                        ..., np.newaxis
                    ],
                )

                # set flow bins to nearest bins
                h_eff = hh.set_flow(h_eff)

                return h_eff

            eff_id = make_efficiency_hist("newveto")
            eff_reco = make_efficiency_hist("reco")
            eff_track = make_efficiency_hist("tracking")

            eff_helper["id"] = makeCorrectionsTensor(eff_id)
            eff_helper["reco"] = makeCorrectionsTensor(eff_reco)
            eff_helper["track"] = makeCorrectionsTensor(eff_track)


hltString = muon_selections.hlt_string(era)
isoBranch = muon_selections.getIsoBranch(args.isolationDefinition)


def build_graph(df, dataset):
    logger.info(f"build graph for dataset: {dataset.name}")
    results = []

    cvh_helper = data_calibration_helper if dataset.is_data else mc_calibration_helper
    jpsi_helper = data_jpsi_crctn_helper if dataset.is_data else mc_jpsi_crctn_helper

    df = df.Define("weight", "std::copysign(1.0, genWeight)")
    # weightsum = df.SumAndCount("weight")
    df = df.Define("isEvenEvent", "event % 2 == 0")
    df = df.Define("unity", "1.0")
    weightsum = df.SumAndCount("unity")

    df = muon_calibration.define_corrected_muons(
        df, cvh_helper, jpsi_helper, args, dataset, smearing_helper, bias_helper
    )

    # gen match to bare muons to select only prompt muons from MC processes, but also including tau decays
    # status flags in NanoAOD: https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv9/2016ULpostVFP/doc_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1.html
    postFSRmuonDef = f"GenPart_status == 1 && (GenPart_statusFlags & 1 || GenPart_statusFlags & (1 << 5)) && abs(GenPart_pdgId) == 13 && GenPart_pt > {args.vetoGenPartPt}"

    # for muon drop method, select events with two postFSR muons
    df = df.Define("postFSRmuons", postFSRmuonDef)
    df = df.Filter(
        "Sum(postFSRmuons) == 2 && (GenPart_pdgId[postFSRmuons][0] + GenPart_pdgId[postFSRmuons][1] == 0)"
    )
    # positive muon
    df = df.Define(
        "postFSRmuon_pt0",
        "GenPart_pdgId[postFSRmuons][0] < 0 ? GenPart_pt[postFSRmuons][0] : GenPart_pt[postFSRmuons][1]",
    )
    df = df.Define(
        "postFSRmuon_eta0",
        "GenPart_pdgId[postFSRmuons][0] < 0 ? GenPart_eta[postFSRmuons][0] : GenPart_eta[postFSRmuons][1]",
    )
    df = df.Define(
        "postFSRmuon_phi0",
        "GenPart_pdgId[postFSRmuons][0] < 0 ? GenPart_phi[postFSRmuons][0] : GenPart_phi[postFSRmuons][1]",
    )
    df = df.Define(
        "postFSRmuon_charge0",  # "-1 * std::copysign(1.0, GenPart_pdgId[postFSRmuons][0])"
        "GenPart_pdgId[postFSRmuons][0] < 0 ? -1 : 1",
    )

    # negative muon
    df = df.Define(
        "postFSRmuon_pt1",
        "GenPart_pdgId[postFSRmuons][0] < 0 ? GenPart_pt[postFSRmuons][1] : GenPart_pt[postFSRmuons][0]",
    )
    df = df.Define(
        "postFSRmuon_eta1",
        "GenPart_pdgId[postFSRmuons][0] < 0 ? GenPart_eta[postFSRmuons][1] : GenPart_eta[postFSRmuons][0]",
    )
    df = df.Define(
        "postFSRmuon_phi1",
        "GenPart_pdgId[postFSRmuons][0] < 0 ? GenPart_phi[postFSRmuons][1] : GenPart_phi[postFSRmuons][0]",
    )
    df = df.Define(
        "postFSRmuon_charge1",  # "-1 * std::copysign(1.0, GenPart_pdgId[postFSRmuons][1])"
        "GenPart_pdgId[postFSRmuons][0] < 0 ? 1 : -1",
    )
    results.append(
        df.HistoBoost(
            f"drop",
            [axis_eta0, axis_pt0, axis_eta1, axis_pt1],
            [
                "postFSRmuon_eta0",
                "postFSRmuon_pt0",
                "postFSRmuon_eta1",
                "postFSRmuon_pt1",
            ],
        )
    )  # [*cols, "nominal_weight"])

    if args.selectReco:
        # define veto muons
        df = muon_selections.select_veto_muons(
            df,
            nMuons=0,
            condition=">=",
            ptCut=15,
            etaCut=2.4,
            useGlobalOrTrackerVeto=False,
            tightGlobalOrTracker=True,
        )
        df = df.Define("oneOrMoreVetoMuons", "Sum(vetoMuons) > 0")

        df = muon_selections.select_good_muons(
            df,
            26,
            56,
            dataset.group,
            nMuons=0,
            condition=">=",
        )
        df = df.Define(
            "GoodTrigObjs",
            f"wrem::goodMuonTriggerCandidate<wrem::Era::Era_{era}>(TrigObj_id,TrigObj_filterBits)",
        )
        df = df.Define(
            "passTrigger",
            f"goodMuons && {hltString} && wrem::hasTriggerMatch(Muon_eta,Muon_phi,TrigObj_eta[GoodTrigObjs],TrigObj_phi[GoodTrigObjs])",
        )
        df = df.Define(
            "passIso",
            f"passTrigger && {isoBranch} < {args.isolationThreshold}",
        )

        # for muon veto efficiency
        for idx in (0, 1):
            odx = 1 if idx == 0 else 0
            # efficiencies measured under the condition that the first muon passes the more stringent cuts
            # df_muon = df.Filter(f"oneOrMoreVetoMuons && wrem::hasMatchDR2(postFSRmuon_eta{idx},postFSRmuon_phi{idx},Muon_eta[vetoMuons],Muon_phi[vetoMuons],0.09)")
            df_muon = df.Filter(
                f"oneOrMoreVetoMuons && wrem::hasMatchDR2(postFSRmuon_eta{idx},postFSRmuon_phi{idx},Muon_eta[passIso],Muon_phi[passIso],0.09)"
            )

            # df_muon = df_muon.Filter(f"postFSRmuon_pt{idx} > 26 && postFSRmuon_pt{idx} < 56")
            # might have more veto muons, but will look for at least one gen matched to the only gen muon
            df_muon = df_muon.Define(
                "passVeto",
                f"oneOrMoreVetoMuons && wrem::hasMatchDR2(postFSRmuon_eta{odx},postFSRmuon_phi{odx},Muon_eta[vetoMuons],Muon_phi[vetoMuons],0.09)",
            )

            results.append(
                df_muon.HistoBoost(
                    f"vetoMuonPass{odx}",
                    [axis_eta0, axis_pt0, axis_passVeto],
                    [f"postFSRmuon_eta{odx}", f"postFSRmuon_pt{odx}", "passVeto"],
                )
            )

    if args.makeCorr:
        return results, weightsum

    df = df.Define(
        "inAcc0",
        f"postFSRmuon_pt0 > {veto_minpt} && abs(postFSRmuon_eta0) < {veto_maxeta}",  # && postFSRmuon_pt0 < {veto_maxpt}
    )
    df = df.Define(
        "inAcc1",
        f"postFSRmuon_pt1 > {veto_minpt} && abs(postFSRmuon_eta1) < {veto_maxeta}",  # && postFSRmuon_pt1 < {veto_maxpt}
    )

    if args.selectReco:
        df = df.Define(
            "muIdx0",
            "wrem::hasMatchDR2idx_closest(postFSRmuon_eta0, postFSRmuon_phi0, Muon_eta[vetoMuons], Muon_phi[vetoMuons], 0.09)",
        )
        df = df.Define(
            "muIdx1",
            "wrem::hasMatchDR2idx_closest(postFSRmuon_eta1, postFSRmuon_phi1, Muon_eta[vetoMuons], Muon_phi[vetoMuons], 0.09)",
        )
        # # make sure also the charge is correct
        df = df.Define(
            "muonIdx0",
            "Muon_charge[vetoMuons][muIdx0] * postFSRmuon_charge0 > 0 ? muIdx0 : -1",
        )
        df = df.Define(
            "muonIdx1",
            "Muon_charge[vetoMuons][muIdx1] * postFSRmuon_charge1 > 0 ? muIdx1 : -1",
        )

    # single muon events (observed)
    df_s = df
    if args.selectReco:
        df_s = df_s.Filter("(muonIdx0 != -1) != (muonIdx1 != -1)")

    def add_hists_obs(df, idx=0, suffix="drop"):
        if args.selectReco:
            df = df.Define(
                f"vtxAgnRelIso04All{idx}",
                f"Muon_vtxAgnPfRelIso04_all[vetoMuons][idx{idx}]",
            )
            df = df.Define(
                f"relIso04All{idx}", f"Muon_pfRelIso04_all[vetoMuons][idx{idx}]"
            )
            df = df.Define(f"dxy{idx}", f"Muon_dxy[vetoMuons][idx{idx}]")
            df = df.Define(f"dxybs{idx}", f"Muon_dxybs[vetoMuons][idx{idx}]")
            results.append(
                df.HistoBoost(
                    f"{suffix}_relIso_vtxAgn04All",
                    [axis_eta, axis_pt, axis_charge, axis_relIso],
                    [
                        f"eta{idx}",
                        f"pt{idx}",
                        f"charge{idx}",
                        f"vtxAgnRelIso04All{idx}",
                    ],
                )
            )
            results.append(
                df.HistoBoost(
                    f"{suffix}_relIso",
                    [axis_eta, axis_pt, axis_charge, axis_relIso],
                    [f"eta{idx}", f"pt{idx}", f"charge{idx}", f"relIso04All{idx}"],
                )
            )
            results.append(
                df.HistoBoost(
                    f"{suffix}_dxy",
                    [axis_eta, axis_pt, axis_charge, axis_dxy],
                    [f"eta{idx}", f"pt{idx}", f"charge{idx}", f"dxy{idx}"],
                )
            )
            results.append(
                df.HistoBoost(
                    f"{suffix}_dxy_bs",
                    [axis_eta, axis_pt, axis_charge, axis_dxy],
                    [f"eta{idx}", f"pt{idx}", f"charge{idx}", f"dxybs{idx}"],
                )
            )

        results.append(
            df.HistoBoost(
                f"{suffix}_met",
                [axis_eta, axis_pt, axis_charge, axis_met],
                [f"eta{idx}", f"pt{idx}", f"charge{idx}", "DeepMETResolutionTune_pt"],
            )
        )
        df = df.Define(
            f"transverseMass{idx}",
            f"wrem::mt_2(pt{idx}, phi{idx}, DeepMETResolutionTune_pt, DeepMETResolutionTune_phi)",
        )
        results.append(
            df.HistoBoost(
                f"{suffix}_mt",
                [axis_eta, axis_pt, axis_charge, axis_mt],
                [f"eta{idx}", f"pt{idx}", f"charge{idx}", f"transverseMass{idx}"],
            )
        )
        df = df.Define(
            f"deltaPhiMuonMet{idx}",
            f"std::abs(wrem::deltaPhi(phi{idx}, DeepMETResolutionTune_phi))",
        )
        results.append(
            df.HistoBoost(
                f"{suffix}_phi",
                [axis_eta, axis_pt, axis_charge, axis_phi],
                [f"eta{idx}", f"pt{idx}", f"charge{idx}", f"deltaPhiMuonMet{idx}"],
            )
        )

    # case A) both muons inside acceptance
    if args.selectReco:
        df_a = df_s.Filter("inAcc0 && inAcc1")
        # check if one muon has a veto muon match and the other one has not
        df_a = df_a.Filter("(muonIdx0 != -1) != (muonIdx1 != -1)")

        df_a = df_a.Define(f"idx0", f"(muonIdx0 != -1) ? muonIdx0 : muonIdx1")
        df_a = df_a.Filter(f"passTrigger[vetoMuons][idx0]")

        df_a = df_a.Define(
            f"pt0", f"(muonIdx0 != -1) ? postFSRmuon_pt0 : postFSRmuon_pt1"
        )
        # df_a = df_a.Filter(f"pt0 > 26 && pt0 < 56")

        df_a = df_a.Define(
            f"eta0", f"(muonIdx0 != -1) ? postFSRmuon_eta0 : postFSRmuon_eta1"
        )
        df_a = df_a.Define(
            f"phi0", f"(muonIdx0 != -1) ? postFSRmuon_phi0 : postFSRmuon_phi1"
        )
        df_a = df_a.Define(
            f"charge0", f"(muonIdx0 != -1) ? postFSRmuon_charge0 : postFSRmuon_charge1"
        )

        add_hists_obs(df_a, suffix="veto")

    # case B) one muon inside accepatance
    df_b = df_s.Filter("inAcc0 != inAcc1")
    if args.selectReco:
        # check if gen muon has a veto muon match and the other one has not
        df_b = df_b.Define("idx0", "inAcc0 ? muonIdx0 : muonIdx1")
        df_b = df_b.Define("idx1", "inAcc0 ? muonIdx1 : muonIdx0")
        df_b = df_b.Filter(f"(idx0 != -1) && (idx1 == -1)")
        df_b = df_b.Filter(f"passTrigger[vetoMuons][idx0]")

    df_b = df_b.Define("pt0", "inAcc0 ? postFSRmuon_pt0 : postFSRmuon_pt1")
    # df_b = df_b.Filter("pt0 > 26 && pt0 < 56")

    df_b = df_b.Define("eta0", "inAcc0 ? postFSRmuon_eta0 : postFSRmuon_eta1")
    df_b = df_b.Define("phi0", "inAcc0 ? postFSRmuon_phi0 : postFSRmuon_phi1")
    df_b = df_b.Define("charge0", "inAcc0 ? postFSRmuon_charge0 : postFSRmuon_charge1")

    add_hists_obs(df_b)

    # select events with both muons inside veto accepatance for closure test
    df_d = df.Filter("inAcc0 && inAcc1")

    if args.selectReco:
        # check if both gen muons have a veto muon match
        df_d = df_d.Filter("(muonIdx0 != -1) && (muonIdx1 != -1)")

        df_d = df_d.Define(
            "vtxAgnRelIso04All0", "Muon_vtxAgnPfRelIso04_all[vetoMuons][muonIdx0]"
        )
        df_d = df_d.Define("relIso04All0", "Muon_pfRelIso04_all[vetoMuons][muonIdx0]")
        df_d = df_d.Define("dxy0", "Muon_dxy[vetoMuons][muonIdx0]")
        df_d = df_d.Define("dxybs0", "Muon_dxybs[vetoMuons][muonIdx0]")

        df_d = df_d.Define(
            "vtxAgnRelIso04All1", "Muon_vtxAgnPfRelIso04_all[vetoMuons][muonIdx1]"
        )
        df_d = df_d.Define("relIso04All1", "Muon_pfRelIso04_all[vetoMuons][muonIdx1]")
        df_d = df_d.Define("dxy1", "Muon_dxy[vetoMuons][muonIdx1]")
        df_d = df_d.Define("dxybs1", "Muon_dxybs[vetoMuons][muonIdx1]")

        if args.tnpEfficiency:
            # for reco efficiencies
            df_d = df_d.Define(
                "trk_pt0", "static_cast<double>(Muon_pt[vetoMuons][muonIdx0])"
            )
            df_d = df_d.Define(
                "trk_pt1", "static_cast<double>(Muon_pt[vetoMuons][muonIdx1])"
            )
            df_d = df_d.Define(
                "trk_eta0", "static_cast<double>(Muon_eta[vetoMuons][muonIdx0])"
            )
            df_d = df_d.Define(
                "trk_eta1", "static_cast<double>(Muon_eta[vetoMuons][muonIdx1])"
            )
            df_d = df_d.Define("trk_charge0", "Muon_charge[vetoMuons][muonIdx0]")
            df_d = df_d.Define("trk_charge1", "Muon_charge[vetoMuons][muonIdx1]")
            # for track efficiency
            df_d = df_d.Define(
                "sa_pt0", "static_cast<double>(Muon_standalonePt[vetoMuons][muonIdx0])"
            )
            df_d = df_d.Define(
                "sa_pt1", "static_cast<double>(Muon_standalonePt[vetoMuons][muonIdx1])"
            )
            df_d = df_d.Define(
                "sa_eta0",
                "static_cast<double>(Muon_standaloneEta[vetoMuons][muonIdx0])",
            )
            df_d = df_d.Define(
                "sa_eta1",
                "static_cast<double>(Muon_standaloneEta[vetoMuons][muonIdx1])",
            )
            df_d = df_d.Define("sa_charge0", "Muon_charge[vetoMuons][muonIdx0]")
            df_d = df_d.Define("sa_charge1", "Muon_charge[vetoMuons][muonIdx1]")

    df_d = df_d.Define("eta0", "static_cast<double>(postFSRmuon_eta0)")
    df_d = df_d.Define("eta1", "static_cast<double>(postFSRmuon_eta1)")
    df_d = df_d.Define("pt0", "static_cast<double>(postFSRmuon_pt0)")
    df_d = df_d.Define("pt1", "static_cast<double>(postFSRmuon_pt1)")
    df_d = df_d.Define("phi0", "static_cast<double>(postFSRmuon_phi0)")
    df_d = df_d.Define("phi1", "static_cast<double>(postFSRmuon_phi1)")
    df_d = df_d.Define(f"charge0", "1")
    df_d = df_d.Define(f"charge1", "-1")

    for idx in (0, 1):
        odx = 1 if idx == 0 else 0

        # select evetns where positive/negative muon is the good one
        df_muon = df_d.Filter(f"passTrigger[vetoMuons][muonIdx{idx}]")
        # df_muon = df_d.Filter(f"postFSRmuon_pt{idx} > 26 && postFSRmuon_pt{idx} < 56")

        if args.selectReco:

            if args.tnpEfficiency:
                df_muon = df_muon.Define(
                    f"eff_reco_weight{idx}_tensor",
                    eff_helper["reco"],
                    [
                        f"trk_charge{odx}",
                        f"trk_eta{odx}",
                        f"trk_pt{odx}",
                        "unity",
                    ],
                )
                df_muon = df_muon.Define(
                    f"eff_id_weight{idx}_tensor",
                    eff_helper["id"],
                    [
                        f"trk_charge{odx}",
                        f"trk_eta{odx}",
                        f"trk_pt{odx}",
                        "unity",
                    ],
                )
                df_muon = df_muon.Define(
                    f"eff_track_weight{idx}_tensor",
                    eff_helper["track"],
                    [
                        f"sa_charge{odx}",
                        f"sa_eta{odx}",
                        f"sa_pt{odx}",
                        "unity",
                    ],
                )
                df_muon = df_muon.Define(
                    f"eff_weight{idx}",
                    f"eff_reco_weight{idx}_tensor[0] * eff_id_weight{idx}_tensor[0] * eff_track_weight{idx}_tensor[0]",
                )

            else:
                df_muon = df_muon.Define(
                    f"eff_weight{idx}_tensor",
                    eff_helper[dataset.name],
                    [
                        f"charge{odx}",
                        f"eta{odx}",
                        f"pt{odx}",
                        "unity",
                    ],
                )
                df_muon = df_muon.Define(
                    f"eff_weight{idx}", f"eff_weight{idx}_tensor[1]"
                )

        # use drop method for event from high mass sample
        df_muon = df_muon.Define(
            f"drop_weight{idx}_tensor",
            drop_helper,
            [
                f"eta{idx}",
                f"pt{idx}",
                f"charge{idx}",
                "unity",
            ],
        )
        if args.selectReco:
            # correct for muon that is dropped
            df_muon = df_muon.Define(
                f"nominal_weight{idx}",
                f"1/eff_weight{idx} * drop_weight{idx}_tensor[0]",
            )
        else:
            df_muon = df_muon.Define(
                f"nominal_weight{idx}", f"drop_weight{idx}_tensor[0]"
            )

        # use anti-veto method for events from low mass sample
        df_muon = df_muon.Define(
            f"anti_veto_weight{idx}", f"(1 - eff_weight{idx})/eff_weight{idx}"
        )

        df_muon = df_muon.Define(
            f"transverseMass{idx}",
            f"wrem::mt_2(pt{idx}, phi{idx}, DeepMETResolutionTune_pt, DeepMETResolutionTune_phi)",
        )
        df_muon = df_muon.Define(
            f"deltaPhiMuonMet{idx}",
            f"std::abs(wrem::deltaPhi(phi{idx}, DeepMETResolutionTune_phi))",
        )

        for key, weight in [
            ("drop", f"nominal_weight{idx}"),
            ("veto", f"anti_veto_weight{idx}"),
        ]:

            results.append(
                df_muon.HistoBoost(
                    f"pred_{key}_met{idx}",
                    [axis_eta, axis_pt, axis_met],
                    [
                        f"eta{idx}",
                        f"pt{idx}",
                        "DeepMETResolutionTune_pt",
                        weight,
                    ],
                )
            )

            results.append(
                df_muon.HistoBoost(
                    f"pred_{key}_mt{idx}",
                    [axis_eta, axis_pt, axis_mt],
                    [
                        f"eta{idx}",
                        f"pt{idx}",
                        f"transverseMass{idx}",
                        weight,
                    ],
                )
            )

            results.append(
                df_muon.HistoBoost(
                    f"pred_{key}_phi{idx}",
                    [axis_eta, axis_pt, axis_phi],
                    [
                        f"eta{idx}",
                        f"pt{idx}",
                        f"deltaPhiMuonMet{idx}",
                        weight,
                    ],
                )
            )

            if args.selectReco:
                results.append(
                    df_muon.HistoBoost(
                        f"pred_{key}_relIso_vtxAgn04All{idx}",
                        [axis_eta, axis_pt, axis_relIso],
                        [
                            f"eta{idx}",
                            f"pt{idx}",
                            f"vtxAgnRelIso04All{idx}",
                            weight,
                        ],
                    )
                )
                results.append(
                    df_muon.HistoBoost(
                        f"pred_{key}_relIso{idx}",
                        [axis_eta, axis_pt, axis_relIso],
                        [
                            f"eta{idx}",
                            f"pt{idx}",
                            f"relIso04All{idx}",
                            weight,
                        ],
                    )
                )
                results.append(
                    df_muon.HistoBoost(
                        f"pred_{key}_dxy{idx}",
                        [axis_eta, axis_pt, axis_dxy],
                        [
                            f"eta{idx}",
                            f"pt{idx}",
                            f"dxy{idx}",
                            weight,
                        ],
                    )
                )
                results.append(
                    df_muon.HistoBoost(
                        f"pred_{key}_dxy_bs{idx}",
                        [axis_eta, axis_pt, axis_dxy],
                        [
                            f"eta{idx}",
                            f"pt{idx}",
                            f"dxybs{idx}",
                            weight,
                        ],
                    )
                )

        # add other muon to MET
        df_muon = df_muon.Define(
            "met_mom4",
            f"ROOT::Math::PtEtaPhiMVector(DeepMETResolutionTune_pt, 0, DeepMETResolutionTune_phi, 0)",
        )
        df_muon = df_muon.Define(
            "muon_mom4",
            f"ROOT::Math::PtEtaPhiMVector(0.1*postFSRmuon_pt{odx}, postFSRmuon_eta{odx}, postFSRmuon_phi{odx}, wrem::muon_mass)",
        )
        df_muon = df_muon.Define(
            f"met{idx}_mom4",
            f"ROOT::Math::PxPyPzEVector(met_mom4)+ROOT::Math::PxPyPzEVector(muon_mom4)",
        )
        df_muon = df_muon.Define(f"met{idx}_pt", f"met{idx}_mom4.pt()")
        df_muon = df_muon.Define(f"met{idx}_phi", f"met{idx}_mom4.phi()")

        results.append(
            df_muon.HistoBoost(
                f"pred_drop_new_met{idx}",
                [axis_eta, axis_pt, axis_met],
                [f"eta{idx}", f"pt{idx}", f"met{idx}_pt", f"nominal_weight{idx}"],
            )
        )
        df_muon = df_muon.Define(
            f"transverseMass{idx}_new",
            f"wrem::mt_2(pt{idx}, phi{idx}, met{idx}_pt, met{idx}_phi)",
        )
        results.append(
            df_muon.HistoBoost(
                f"pred_drop_new_mt{idx}",
                [axis_eta, axis_pt, axis_mt],
                [
                    f"eta{idx}",
                    f"pt{idx}",
                    f"transverseMass{idx}_new",
                    f"nominal_weight{idx}",
                ],
            )
        )
        df_muon = df_muon.Define(
            f"deltaPhiMuonMet{idx}_new",
            f"std::abs(wrem::deltaPhi(phi{idx}, met{idx}_phi))",
        )
        results.append(
            df_muon.HistoBoost(
                f"pred_drop_new_phi{idx}",
                [axis_eta, axis_pt, axis_phi],
                [
                    f"eta{idx}",
                    f"pt{idx}",
                    f"deltaPhiMuonMet{idx}_new",
                    f"nominal_weight{idx}",
                ],
            )
        )

    # # for muon veto transfer method
    # for sign in ("plus", "minus"):
    #     # select
    #     df_muon = df.Define(
    #         "postFSRmuons",
    #         f"{postFSRmuonDef} && GenPart_pdgId=={'13' if sign == 'plus' else '-13'}",
    #     )
    #     # df = df.Define("postFSRantimuons", f"{postFSRmuonDef} && GenPart_pdgId==-13")

    #     # restrict to one gen muon for simplicity
    #     df_muon = df_muon.Filter("Sum(postFSRmuons) == 1")
    #     # df = muon_selections.veto_electrons(df)
    #     # df = muon_selections.apply_met_filters(df)

    #     ########################################################################
    #     # # define event weights here since they are needed below for some helpers
    #     # df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
    #     # df = df.Define("weight_vtx", vertex_helper, ["GenVtx_z", "Pileup_nTrueInt"])

    #     # weight_expr = "weight_pu*L1PreFiringWeight_ECAL_Nom"
    #     # if not args.noVertexWeight:
    #     #     weight_expr += "*weight_vtx"

    #     # logger.debug(f"Exp weight defined: {weight_expr}")
    #     # df = df.Define("exp_weight", weight_expr)
    #     # df = theory_tools.define_theory_weights_and_corrs(
    #     #     df, dataset.name, corr_helpers, args
    #     # )

    #     ########################################################################

    #     df_muon = df_muon.Define("postFSRmuon_pt0", "GenPart_pt[postFSRmuons][0]")
    #     df_muon = df_muon.Define("postFSRmuon_eta0", "GenPart_eta[postFSRmuons][0]")
    #     df_muon = df_muon.Define("postFSRmuon_phi0", "GenPart_phi[postFSRmuons][0]")
    #     df_muon = df_muon.Define(
    #         "postFSRmuon_charge0",
    #         "-1 * std::copysign(1.0, GenPart_pdgId[postFSRmuons][0])",
    #     )

    #     df_muon = muon_selections.select_veto_muons(
    #         df_muon,
    #         nMuons=0,
    #         condition=">=",
    #         ptCut=10,  # args.vetoRecoPt,
    #         etaCut=3.0,
    #         useGlobalOrTrackerVeto=False,
    #         tightGlobalOrTracker=True,
    #     )
    #     # might have more veto muons, but will look for at least one gen matched to the only gen muon
    #     df_muon = df_muon.Define("oneOrMoreVetoMuons", "Sum(vetoMuons) > 0")
    #     df_muon = df_muon.Define(
    #         "passVeto",
    #         "oneOrMoreVetoMuons && wrem::hasMatchDR2(postFSRmuon_eta0,postFSRmuon_phi0,Muon_eta[vetoMuons],Muon_phi[vetoMuons],0.09)",
    #     )

    #     nominal = df_muon.HistoBoost(
    #         f"nominal_{sign}", axes, cols
    #     )  # [*cols, "nominal_weight"])
    #     results.append(nominal)

    return results, weightsum


resultdict = narf.build_and_run(datasets, build_graph)

if not args.noScaleToData:
    scale_to_data(resultdict)
    aggregate_groups(datasets, resultdict, groups_to_aggregate)

write_analysis_output(
    resultdict, f"{os.path.basename(__file__).replace('py', 'hdf5')}", args
)
