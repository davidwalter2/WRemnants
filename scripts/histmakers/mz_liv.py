import argparse
import os
import narf
import hist

import pdb

from utilities import logging, common
from utilities.io_tools import output_tools
from wremnants import pileup, vertex, muon_prefiring, syst_tools
from wremnants.datasets.dataset_tools import getDatasets

from zcount import muon_utils
from zcount.lumi_tools import get_run_lumi_axes

parser,initargs = common.common_parser(True)
args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose)

# definitions
ptLow = 25
ptHigh = 200

ptLowStandalone = 25
ptHighStandalone = 200

massGenMin = 66
massGenMax = 116

nStandaloneValidHist = 1

axis_mll = hist.axis.Regular(60, 60., 120., name = "mll", overflow=True, underflow=True)
axis_number = hist.axis.Integer(0, 10, name = "number", overflow=True, underflow=False)
axis_nPU = hist.axis.Integer(0, 100, name = "nPU", overflow=True, underflow=False)
axis_nPV = hist.axis.Integer(0, 100, name = "nPV", overflow=True, underflow=False)
axis_acceptance_gen = hist.axis.Integer(-1, 4, name = "acceptance_gen", overflow=False, underflow=False)
axis_acceptance_reco = hist.axis.Integer(0, 4, name = "acceptance_reco", overflow=False, underflow=False)
axis_os = hist.axis.Boolean(name = "os")

column_nPU = "Pileup_nTrueInt"
column_nPV = "PV_npvsGood"

axes_n = [axis_number]

# get datasets
datasets = getDatasets(
    base_path=args.dataPath,
    maxFiles=args.maxFiles,
    filt=args.filterProcs,
    excl=args.excludeProcs)

tnpNano = False

# initialize helpers
# pileup_helper = pileup.make_pileup_helper(era = args.era)
vertex_helper = vertex.make_vertex_helper(era = args.era)
muon_prefiring_helper, muon_prefiring_helper_stat, muon_prefiring_helper_syst = muon_prefiring.make_muon_prefiring_helpers(era = args.era)

# helper function
def add_hist_with_syst(df, results, name, axes, cols, dataset):
    if dataset.is_data:
        logger.info("Add data hist")
        results.append(df.HistoBoost(name, axes, cols))
    else:
        logger.info("Add MC hist")
        results.append(df.HistoBoost(name, axes, [*cols, "nominal_weight"]))
        syst_tools.add_L1Prefire_unc_hists(results, df, muon_prefiring_helper_stat, muon_prefiring_helper_syst, axes, cols, base_name=name)


def build_graph(df, dataset):
    logger.info(f"build graph for dataset: {dataset.name}")
    results = []

    isZ = dataset.name in ["ZmumuPostVFP",]

    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    weightsum = df.SumAndCount("weight")

    ### TODO: muon calibration

    ### muon selection
    if not isZ:
        df = df.Filter("(HLT_IsoTkMu24 || HLT_IsoMu24)")        
    else:
        # gen level matching
        # use post FSR muons, including stable muons from tau decays
        df = df.Define("postFSRmuons", "GenPart_status == 1 && (GenPart_statusFlags & 1 || GenPart_statusFlags & (5<<1)) && GenPart_pdgId == 13")
        df = df.Define("postFSRantimuons", "GenPart_status == 1 && (GenPart_statusFlags & 1 || GenPart_statusFlags & (5<<1)) && GenPart_pdgId == -13")
        # in case of multiple gen muons, take the ones with the highest pT
        df = df.Define("postFSRmuonIdx", "ROOT::VecOps::ArgMax(GenPart_pt[postFSRmuons])")
        df = df.Define("postFSRantimuonIdx", "ROOT::VecOps::ArgMax(GenPart_pt[postFSRantimuons])")

        df = df.Define("GenMuon_pt", f"GenPart_pt[postFSRmuons][postFSRmuonIdx]")
        df = df.Define("GenAntiMuon_pt", f"GenPart_pt[postFSRantimuons][postFSRantimuonIdx]")
        df = df.Define("GenMuon_eta", f"GenPart_eta[postFSRmuons][postFSRmuonIdx]")
        df = df.Define("GenAntiMuon_eta", f"GenPart_eta[postFSRantimuons][postFSRantimuonIdx]")
        # define some dummy columns to run code more harmonized
        df = df.DefinePerSample("GenMuon_charge", f"-1")
        df = df.DefinePerSample("GenAntiMuon_charge", f"1")
        df = df.DefinePerSample("GenMuon_isGen", f"1")
        df = df.DefinePerSample("GenAntiMuon_isGen", f"1")

        df = df.Define("GenMuon_mom4", "ROOT::Math::PtEtaPhiMVector(GenPart_pt[postFSRmuons][postFSRmuonIdx], GenPart_eta[postFSRmuons][postFSRmuonIdx], GenPart_phi[postFSRmuons][postFSRmuonIdx], GenPart_mass[postFSRmuons][postFSRmuonIdx])")
        df = df.Define("GenAntiMuon_mom4", "ROOT::Math::PtEtaPhiMVector(GenPart_pt[postFSRantimuons][postFSRantimuonIdx], GenPart_eta[postFSRantimuons][postFSRantimuonIdx], GenPart_phi[postFSRantimuons][postFSRantimuonIdx], GenPart_mass[postFSRantimuons][postFSRantimuonIdx])")
        df = df.Define("GenZ_mom4", "ROOT::Math::PxPyPzEVector(GenMuon_mom4)+ROOT::Math::PxPyPzEVector(GenAntiMuon_mom4)")
        df = df.Define("GenZ_mass", "GenZ_mom4.mass()")
 
        # gen level acceptance on an event level (0: out of acceptance; 1: endcap-endcap; 2: endcap-barrel; 3: barrel-barrel)
        df = df.Define("Event_acceptance_gen", f"""
            (GenMuon_pt > {ptLow} && GenAntiMuon_pt > {ptLow} && GenMuon_pt < {ptHigh} && GenAntiMuon_pt < {ptHigh}
            && std::fabs(GenMuon_eta) < 2.4 && std::fabs(GenAntiMuon_eta) < 2.4
            && GenZ_mass > {massGenMin} && GenZ_mass < {massGenMax})*(1 + (std::fabs(GenMuon_eta) < 0.9) + (std::fabs(GenAntiMuon_eta) < 0.9))
            """)

    # reco acceptance definition (0: out of acceptance; 0.5: endcap; 1.5: barrel)
    df = df.Define("Muon_acceptance", f"(Muon_pt > {ptLow} && Muon_pt < {ptHigh} && abs(Muon_eta) < 2.4)*(0.5 + (abs(Muon_eta) < 0.9))")

    # minimal selection for protection (only objects with pT>15GeV are stored independently)
    baseline_muon_selection = "&& Muon_pt > 15"# && abs(Muon_eta) < 2.4"
    baseline_muon_standalone_selection = "&& Muon_standalonePt > 15"# && abs(Muon_standaloneEta) < 2.4"

    if not isZ:
        acceptance_selection = "&& Muon_acceptance > 0"
    else:
        acceptance_selection = baseline_muon_selection

    if tnpNano:
        # info only available in tnp nanoAOD
        df = df.Define("Track_acceptance", f"(Track_pt > {ptLow} && Track_pt < {ptHigh} && abs(Track_eta) < 2.4)*(0.5 + (abs(Track_eta) < 0.9))")
        df = df.Define("MergedStandAloneMuon_acceptance", f"""
            (MergedStandAloneMuon_pt > {ptLowStandalone} && MergedStandAloneMuon_pt < {ptHighStandalone} 
            && abs(MergedStandAloneMuon_eta) < 2.4)*(0.5 + (abs(MergedStandAloneMuon_eta) < 0.9))""")

        baseline_track_selection = "&& Track_pt > 15"# && abs(Track_eta) < 2.4"
        baseline_standalone_selection = "&& MergedStandAloneMuon_pt > 15"# && abs(MergedStandAloneMuon_eta) < 2.4"

        if not isZ:
            acceptance_track_selection = "&& Track_acceptance > 0"
            acceptance_standalone_selection = "&& MergedStandAloneMuon_acceptance > 0"
        else:
            acceptance_track_selection = baseline_track_selection
            acceptance_standalone_selection = baseline_standalone_selection

        ### track muon definition starting from track collection
        df = df.Define("Track_isGoodTrack", f"Track_trackOriginalAlgo != 13 && Track_trackOriginalAlgo != 14 && (Track_qualityMask & 4) {acceptance_track_selection}")
        
        # matched tracks with collection of good standalone muons
        df = df.Define("MergedStandAloneMuon_isGoodForTrack", f"MergedStandAloneMuon_numberOfValidHits >= {nStandaloneValidHist} {baseline_standalone_selection}")
        df = df.Define("Track_hasGoodStandalone", f"""Track_isGoodTrack
            && zcount::hasAnyMatchDR2(Track_eta, Track_phi, 
                MergedStandAloneMuon_eta[MergedStandAloneMuon_isGoodForTrack], MergedStandAloneMuon_phi[MergedStandAloneMuon_isGoodForTrack], 0.09)
            """)

        #### standalone muon definition starting from outer track
        df = df.Define("MergedStandAloneMuon_isGoodStandalone", f"MergedStandAloneMuon_numberOfValidHits >= {nStandaloneValidHist} {acceptance_standalone_selection}")
        df = df.Define("Muon_isGoodForStandalone", f"Muon_isGlobal && Muon_highPurity {baseline_muon_selection}")
        df = df.Define("MergedStandAloneMuon_isGoodGlobal", f"""MergedStandAloneMuon_isGoodStandalone
            && zcount::mergedStandAloneIsGoodGlobal(
                MergedStandAloneMuon_extraIdx, MergedStandAloneMuon_eta, MergedStandAloneMuon_phi, 
                Muon_standaloneExtraIdx[Muon_isGoodForStandalone], Muon_eta[Muon_isGoodForStandalone], Muon_phi[Muon_isGoodForStandalone])""")

        # some auxiliary histograms
        df = df.Define("nStandalone", "Sum(MergedStandAloneMuon_isGoodStandalone)")
        results.append(df.HistoBoost(f"nStandalone", [axis_number,], ["nStandalone"]))

        # make mutually exclsive categories of four vectors
        df = muon_utils.make_collection(df, name="TrkFailMuon", input_name="Track", selection="Track_isGoodTrack && !Track_hasGoodStandalone", do_gen_match=isZ)
        df = muon_utils.make_collection(df, name="TrkPassMuon", input_name="Track", selection="Track_hasGoodStandalone", do_gen_match=isZ)

        df = muon_utils.make_collection(df, name="StaFailMuon", input_name="MergedStandAloneMuon", selection="MergedStandAloneMuon_isGoodStandalone && !MergedStandAloneMuon_isGoodGlobal", do_gen_match=isZ)
        df = muon_utils.make_collection(df, name="StaPassMuon", input_name="MergedStandAloneMuon", selection="MergedStandAloneMuon_isGoodGlobal", do_gen_match=isZ)

        categories = [
            ("HLT_TrkFail", "HLT", "TrkFail"),
            ("HLT_TrkPass", "HLT", "TrkPass"),
            ("HLT_StaFail", "HLT", "StaFail"),
            ("HLT_StaPass", "HLT", "StaPass"),
            ]
    else:
        categories = []

    #### global muon definition starting from inner track
    if tnpNano:
        muon_nStandaloneValidHits_str = f"zcount::hasTrackRef(Muon_standaloneExtraIdx, MergedStandAloneMuon_extraIdx[MergedStandAloneMuon_numberOfValidHits >= {nStandaloneValidHist}])"
    else:
        muon_nStandaloneValidHits_str = f"Muon_standaloneNumberOfValidHits >= {nStandaloneValidHist}"

    df = df.Define("Muon_isGoodGlobal", f"""Muon_isGlobal && Muon_highPurity {acceptance_selection} {baseline_muon_standalone_selection}
        && {muon_nStandaloneValidHits_str} && zcount::hasMatchDR2(Muon_standaloneEta, Muon_standalonePhi, Muon_eta, Muon_phi, 0.09)""")

    #### ID muon definition
    df = df.Define("Muon_isGoodID", "Muon_isGoodGlobal && Muon_mediumId && abs(Muon_dxybs) < 0.05 && Muon_vtxAgnPfRelIso03_chg < 0.15") 

    #### HLT muon definition for good ID 
    # df = df.Define("goodTrigObjs", "wrem::goodMuonTriggerCandidate(TrigObj_id,TrigObj_pt,TrigObj_l1pt,TrigObj_l2pt,TrigObj_filterBits)")
    df = df.Define("goodTrigObjs", "wrem::goodMuonTriggerCandidate(TrigObj_id,TrigObj_filterBits)")
    df = df.Define("Muon_isGoodHLT", "zcount::hasAnyMatchDR2(Muon_eta[Muon_isGoodID],Muon_phi[Muon_isGoodID],TrigObj_eta[goodTrigObjs],TrigObj_phi[goodTrigObjs], 0.09)") 

    # define event weights
    if not dataset.is_data:
        # define nominal weight
        # df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
        df = df.Define("weight_vtx", vertex_helper, ["GenVtx_z", "Pileup_nTrueInt"])
        #TODO: use corrected
        # df = df.Define("weight_newMuonPrefiringSF", muon_prefiring_helper,["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId"])
        df = df.Define("weight_newMuonPrefiringSF", muon_prefiring_helper,["Muon_eta", "Muon_pt", "Muon_phi", "Muon_charge", "Muon_looseId"])

        weight_expr = "weight_vtx*weight_newMuonPrefiringSF*L1PreFiringWeight_ECAL_Nom"
        df = df.Define("nominal_weight", weight_expr)

    df = muon_utils.make_collection(df, name="GloMuon", selection="Muon_isGoodGlobal && !Muon_isGoodID", do_gen_match=isZ)
    df = muon_utils.make_collection(df, name="IDMuon", selection="Muon_isGoodID", subselection="!Muon_isGoodHLT", do_gen_match=isZ)
    df = muon_utils.make_collection(df, name="HLTMuon", selection="Muon_isGoodID", subselection="Muon_isGoodHLT", do_gen_match=isZ)

    if dataset.is_data:
        # df = df.Define("runNumber", "static_cast<std::string>(run)")
        axis_run, axis_lumi = get_run_lumi_axes(dataset.lumi_csv, dataset.lumi_json)

    categories += [
        ("HLT_Glo", "HLT", "Glo"),
        ("HLT_ID", "HLT", "ID"),
        ("HLT_HLT", "HLT", None)
    ]
    if isZ:
        categories += [
            ("Gen", "Gen", "GenAnti"),
            # ("HLT_Gen", "HLT", "Gen"),
            # ("HLT_GenAnti", "HLT", "GenAnti"),
            ("Glo_Glo", "Glo", None),
            ("ID_Glo", "ID", "Glo"),
            ("ID_ID", "ID", None),            
        ]

    df = df.Alias("Muon_correctedEta", "Muon_eta")
    df = df.Alias("Muon_correctedPt", "Muon_pt")
    df = df.Alias("Muon_correctedPhi", "Muon_phi")
    df = df.Alias("Muon_correctedCharge", "Muon_charge")
    # df = df.Alias("Muon_looseId", 

    # define histograms in mutually exclsive categories of pairs of muons
    for suffix, object1, object2 in categories:
        logger.info(f"Set histograms for category {suffix}")

        add_acceptance_gen = isZ
        add_acceptance_reco = isZ and suffix not in ["Gen", "HLT_Gen", "HLT_GenAnti"]
        add_os_charge = "Sta" not in suffix
        add_mass = True

        df = muon_utils.make_pairs(df, f"pair_{suffix}", f"{object1}Muon", f"{object2}Muon" if object2 else None, 
            define_acceptance_gen=add_acceptance_gen, define_acceptance_reco=add_acceptance_reco, define_os_charge=add_os_charge, define_mass=add_mass)

        # define column and axes for main histograms
        cols = []
        axes = []
        if add_mass:
            cols.append(f"pair_{suffix}_mass")
            axes.append(axis_mll)
        if add_os_charge:
            cols.append(f"pair_{suffix}_os")
            axes.append(axis_os)
        if add_acceptance_gen:
            cols.append(f"pair_{suffix}_acceptance_gen")
            axes.append(axis_acceptance_gen)
        if add_acceptance_reco:
            cols.append(f"pair_{suffix}_acceptance_reco")
            axes.append(axis_acceptance_reco)
        if dataset.is_data:
            axes.append(axis_run)
            cols.append("run")
            axes.append(axis_lumi)
            cols.append("luminosityBlock")

        # define main histograms
        if not dataset.is_data:
            # pileup histogram
            add_hist_with_syst(df, results, f"pair_{suffix}_nPU_mass", [axis_nPU, *axes], [column_nPU, *cols], dataset)            

        # primary vertex histogram
        add_hist_with_syst(df, results, f"pair_{suffix}_nPV_mass", [axis_nPV, *axes], [column_nPV, *cols], dataset)            

    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)
output_tools.write_analysis_output(resultdict, f"{os.path.basename(__file__).replace('py', 'hdf5')}", args)