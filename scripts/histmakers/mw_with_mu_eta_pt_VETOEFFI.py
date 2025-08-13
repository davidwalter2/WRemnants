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
            "pt1": slice(0, hist.overflow, hist.sum),
        }
    ]
    # observerd positive muon inside, other ourside
    hp_out = hh.addHists(hp_all, hp_in, scale2=-1)

    hm = h1[{"pt1": slice(26j, 56j), "eta1": slice(-2.4j, 2.4j)}]
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
    h2 = loadHist("ZmumuPostVFP", "drop")

    h2p = h2[{"pt": slice(26j, 56j), "eta": slice(-2.4j, 2.4j)}]
    h2m = h2[{"pt1": slice(26j, 56j), "eta1": slice(-2.4j, 2.4j)}]
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


def build_graph(df, dataset):
    logger.info(f"build graph for dataset: {dataset.name}")
    results = []
    isW = dataset.name in common.wprocs
    isWmunu = dataset.name in ["WplusmunuPostVFP", "WminusmunuPostVFP"]
    isZ = dataset.name in common.zprocs
    isWorZ = isW or isZ
    isTop = dataset.group == "Top"

    cvh_helper = data_calibration_helper if dataset.is_data else mc_calibration_helper
    jpsi_helper = data_jpsi_crctn_helper if dataset.is_data else mc_jpsi_crctn_helper

    df = df.Define("weight", "std::copysign(1.0, genWeight)")
    # weightsum = df.SumAndCount("weight")
    df = df.Define("isEvenEvent", "event % 2 == 0")
    df = df.Define("unity", "1.0")
    weightsum = df.SumAndCount("unity")

    axes = nominal_axes
    cols = nominal_cols

    df = muon_calibration.define_corrected_muons(
        df, cvh_helper, jpsi_helper, args, dataset, smearing_helper, bias_helper
    )

    # gen match to bare muons to select only prompt muons from MC processes, but also including tau decays
    # status flags in NanoAOD: https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv9/2016ULpostVFP/doc_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1.html
    postFSRmuonDef = f"GenPart_status == 1 && (GenPart_statusFlags & 1 || GenPart_statusFlags & (1 << 5)) && abs(GenPart_pdgId) == 13 && GenPart_pt > {args.vetoGenPartPt}"

    # for muon drop method, select events with two postFSR muons
    df_dimuon = df.Define("postFSRmuons", postFSRmuonDef)
    df_dimuon = df_dimuon.Filter(
        "Sum(postFSRmuons) == 2 && (GenPart_pdgId[postFSRmuons][0] + GenPart_pdgId[postFSRmuons][1] == 0)"
    )
    # positive muon
    df_dimuon = df_dimuon.Define(
        "postFSRmuon_pt0",
        "GenPart_pdgId[postFSRmuons][0] < 0 ? GenPart_pt[postFSRmuons][0] : GenPart_pt[postFSRmuons][1]",
    )
    df_dimuon = df_dimuon.Define(
        "postFSRmuon_eta0",
        "GenPart_pdgId[postFSRmuons][0] < 0 ? GenPart_eta[postFSRmuons][0] : GenPart_eta[postFSRmuons][1]",
    )
    df_dimuon = df_dimuon.Define(
        "postFSRmuon_phi0",
        "GenPart_pdgId[postFSRmuons][0] < 0 ? GenPart_phi[postFSRmuons][0] : GenPart_phi[postFSRmuons][1]",
    )
    df_dimuon = df_dimuon.Define(
        "postFSRmuon_charge0", "-1 * std::copysign(1.0, GenPart_pdgId[postFSRmuons][0])"
    )

    # negative muon
    df_dimuon = df_dimuon.Define(
        "postFSRmuon_pt1",
        "GenPart_pdgId[postFSRmuons][0] < 0 ? GenPart_pt[postFSRmuons][1] : GenPart_pt[postFSRmuons][0]",
    )
    df_dimuon = df_dimuon.Define(
        "postFSRmuon_eta1",
        "GenPart_pdgId[postFSRmuons][0] < 0 ? GenPart_eta[postFSRmuons][1] : GenPart_eta[postFSRmuons][0]",
    )
    df_dimuon = df_dimuon.Define(
        "postFSRmuon_phi1",
        "GenPart_pdgId[postFSRmuons][0] < 0 ? GenPart_phi[postFSRmuons][1] : GenPart_phi[postFSRmuons][0]",
    )
    df_dimuon = df_dimuon.Define(
        "postFSRmuon_charge1", "-1 * std::copysign(1.0, GenPart_pdgId[postFSRmuons][1])"
    )
    results.append(
        df_dimuon.HistoBoost(
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

    # define veto muons
    df_dimuon = muon_selections.select_veto_muons(
        df_dimuon,
        nMuons=0,
        condition=">=",
        ptCut=15,
        etaCut=2.4,
        useGlobalOrTrackerVeto=False,
        tightGlobalOrTracker=True,
    )
    df_dimuon = df_dimuon.Define("oneOrMoreVetoMuons", "Sum(vetoMuons) > 0")

    # for muon veto efficiency
    for idx in (0, 1):
        odx = 1 if idx == 0 else 0
        # efficiencies measured under the condition that the first muon passes the more stringent cuts
        df_muon = df_dimuon.Filter(
            f"postFSRmuon_pt{idx} > 26 && postFSRmuon_pt{idx} < 56"
        )
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

    df_dimuon = df_dimuon.Define(
        "inAcc0",
        f"postFSRmuon_pt0 > {veto_minpt} && postFSRmuon_pt0 < {veto_maxpt} && abs(postFSRmuon_eta0) < {veto_maxeta}",
    )
    df_dimuon = df_dimuon.Define(
        "inAcc1",
        f"postFSRmuon_pt1 > {veto_minpt} && postFSRmuon_pt1 < {veto_maxpt} && abs(postFSRmuon_eta1) < {veto_maxeta}",
    )

    if dataset.name == "DYJetsToMuMuMass10to50PostVFP":
        # single muon events (observed)
        df_s = df_dimuon.Filter("inAcc0 != inAcc1")
        df_s = df_s.Define("eta0", "inAcc0 ? postFSRmuon_eta0 : postFSRmuon_eta1")
        df_s = df_s.Define("pt0", "inAcc0 ? postFSRmuon_pt0 : postFSRmuon_pt1")

        df_s = df_s.Filter("pt0 > 26 && pt0 < 56")

        df_s = df_s.Define("phi0", "inAcc0 ? postFSRmuon_phi0 : postFSRmuon_phi1")
        df_s = df_s.Define(
            "charge0", "inAcc0 ? postFSRmuon_charge0 : postFSRmuon_charge1"
        )

        if args.selectReco:
            # check if gen muon has a veto muon match
            df_s = df_s.Define(
                "muonIdx",
                "wrem::hasMatchDR2idx(eta0,phi0,Muon_eta[vetoMuons],Muon_phi[vetoMuons],0.09)",
            )
            df_s = df_s.Filter("muonIdx != -1")
            df_s = df_s.Define(
                "vtxAgnRelIso04All0", "Muon_vtxAgnPfRelIso04_all[muonIdx]"
            )
            df_s = df_s.Define("relIso04All0", "Muon_pfRelIso04_all[muonIdx]")
            df_s = df_s.Define("dxy0", "Muon_dxy[muonIdx]")
            df_s = df_s.Define("dxybs0", "Muon_dxybs[muonIdx]")
            results.append(
                df_s.HistoBoost(
                    f"drop_relIso_vtxAgn04All",
                    [axis_eta, axis_pt, axis_charge, axis_relIso],
                    ["eta0", "pt0", "charge0", "vtxAgnRelIso04All0"],
                )
            )
            results.append(
                df_s.HistoBoost(
                    f"drop_relIso",
                    [axis_eta, axis_pt, axis_charge, axis_relIso],
                    ["eta0", "pt0", "charge0", "relIso04All0"],
                )
            )
            results.append(
                df_s.HistoBoost(
                    f"drop_dxy",
                    [axis_eta, axis_pt, axis_charge, axis_dxy],
                    ["eta0", "pt0", "charge0", "dxy0"],
                )
            )
            results.append(
                df_s.HistoBoost(
                    f"drop_dxy_bs",
                    [axis_eta, axis_pt, axis_charge, axis_dxy],
                    ["eta0", "pt0", "charge0", "dxybs0"],
                )
            )

        results.append(
            df_s.HistoBoost(
                f"drop_met",
                [axis_eta, axis_pt, axis_charge, axis_met],
                ["eta0", "pt0", "charge0", "DeepMETResolutionTune_pt"],
            )
        )
        df_s = df_s.Define(
            "transverseMass0",
            "wrem::mt_2(pt0, phi0, DeepMETResolutionTune_pt, DeepMETResolutionTune_phi)",
        )
        results.append(
            df_s.HistoBoost(
                f"drop_mt",
                [axis_eta, axis_pt, axis_charge, axis_mt],
                ["eta0", "pt0", "charge0", "transverseMass0"],
            )
        )
        df_s = df_s.Define(
            "deltaPhiMuonMet0",
            "std::abs(wrem::deltaPhi(phi0, DeepMETResolutionTune_phi))",
        )
        results.append(
            df_s.HistoBoost(
                f"drop_phi",
                [axis_eta, axis_pt, axis_charge, axis_phi],
                ["eta0", "pt0", "charge0", "deltaPhiMuonMet0"],
            )
        )

    if dataset.name == "ZmumuPostVFP":
        # select events with both muons inside veto accepatance for closure test
        df_d = df_dimuon.Filter("inAcc0 && inAcc1")

        if args.selectReco:
            # check if gen muons have a veto muon match
            df_d = df_d.Define(
                "muonIdx0",
                "wrem::hasMatchDR2idx(postFSRmuon_eta0, postFSRmuon_phi0, Muon_eta[vetoMuons], Muon_phi[vetoMuons], 0.09)",
            )
            df_d = df_d.Filter("muonIdx0 != -1")
            df_d = df_d.Define(
                "vtxAgnRelIso04All0", "Muon_vtxAgnPfRelIso04_all[muonIdx0]"
            )
            df_d = df_d.Define("relIso04All0", "Muon_pfRelIso04_all[muonIdx0]")
            df_d = df_d.Define("dxy0", "Muon_dxy[muonIdx0]")
            df_d = df_d.Define("dxybs0", "Muon_dxybs[muonIdx0]")

            df_d = df_d.Define(
                "muonIdx1",
                "wrem::hasMatchDR2idx(postFSRmuon_eta1, postFSRmuon_phi1, Muon_eta[vetoMuons], Muon_phi[vetoMuons], 0.09)",
            )
            df_d = df_d.Filter("muonIdx1 != -1")
            df_d = df_d.Define(
                "vtxAgnRelIso04All1", "Muon_vtxAgnPfRelIso04_all[muonIdx1]"
            )
            df_d = df_d.Define("relIso04All1", "Muon_pfRelIso04_all[muonIdx1]")
            df_d = df_d.Define("dxy1", "Muon_dxy[muonIdx1]")
            df_d = df_d.Define("dxybs1", "Muon_dxybs[muonIdx1]")

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

            # select evetns where positive/negative muon is in final accepatnce
            df_muon = df_d.Filter(
                f"postFSRmuon_pt{idx} > 26 && postFSRmuon_pt{idx} < 56"
            )

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
                    f"nominal_weight{idx}",
                    f"1/eff_weight{idx}_tensor[1] * drop_weight{idx}_tensor[0]",
                )
            else:
                df_muon = df_muon.Define(
                    f"nominal_weight{idx}", f"drop_weight{idx}_tensor[0]"
                )

            results.append(
                df_muon.HistoBoost(
                    f"pred_drop_met{idx}",
                    [axis_eta, axis_pt, axis_met],
                    [
                        f"eta{idx}",
                        f"pt{idx}",
                        "DeepMETResolutionTune_pt",
                        f"nominal_weight{idx}",
                    ],
                )
            )
            df_muon = df_muon.Define(
                f"transverseMass{idx}",
                f"wrem::mt_2(pt{idx}, phi{idx}, DeepMETResolutionTune_pt, DeepMETResolutionTune_phi)",
            )
            results.append(
                df_muon.HistoBoost(
                    f"pred_drop_mt{idx}",
                    [axis_eta, axis_pt, axis_mt],
                    [
                        f"eta{idx}",
                        f"pt{idx}",
                        f"transverseMass{idx}",
                        f"nominal_weight{idx}",
                    ],
                )
            )
            df_muon = df_muon.Define(
                f"deltaPhiMuonMet{idx}",
                f"std::abs(wrem::deltaPhi(phi{idx}, DeepMETResolutionTune_phi))",
            )
            results.append(
                df_muon.HistoBoost(
                    f"pred_drop_phi{idx}",
                    [axis_eta, axis_pt, axis_phi],
                    [
                        f"eta{idx}",
                        f"pt{idx}",
                        f"deltaPhiMuonMet{idx}",
                        f"nominal_weight{idx}",
                    ],
                )
            )

            if args.selectReco:
                results.append(
                    df_muon.HistoBoost(
                        f"pred_drop_relIso_vtxAgn04All{idx}",
                        [axis_eta, axis_pt, axis_relIso],
                        [
                            f"eta{idx}",
                            f"pt{idx}",
                            f"vtxAgnRelIso04All{idx}",
                            f"nominal_weight{idx}",
                        ],
                    )
                )
                results.append(
                    df_muon.HistoBoost(
                        f"pred_drop_relIso{idx}",
                        [axis_eta, axis_pt, axis_relIso],
                        [
                            f"eta{idx}",
                            f"pt{idx}",
                            f"relIso04All{idx}",
                            f"nominal_weight{idx}",
                        ],
                    )
                )
                results.append(
                    df_muon.HistoBoost(
                        f"pred_drop_dxy{idx}",
                        [axis_eta, axis_pt, axis_dxy],
                        [
                            f"eta{idx}",
                            f"pt{idx}",
                            f"dxy{idx}",
                            f"nominal_weight{idx}",
                        ],
                    )
                )
                results.append(
                    df_muon.HistoBoost(
                        f"pred_drop_dxy_bs{idx}",
                        [axis_eta, axis_pt, axis_dxy],
                        [
                            f"eta{idx}",
                            f"pt{idx}",
                            f"dxybs{idx}",
                            f"nominal_weight{idx}",
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
