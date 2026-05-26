import os

from wremnants.utilities import binning, common, parsing, samples, theory_utils

analysis_label = common.analysis_label(os.path.basename(__file__))
parser, initargs = parsing.common_parser(analysis_label)

import math

import hist
import numpy as np
import ROOT

import narf
from wremnants.production import (
    muon_calibration,
    muon_efficiencies_binned,
    muon_efficiencies_insitu,
    muon_efficiencies_smooth,
    muon_prefiring,
    muon_selections,
    pileup,
    systematics,
    theory_corrections,
    theoryAgnostic_tools,
    unfolding_tools,
    vertex,
)
from wremnants.production.datasets.dataset_tools import getDatasets
from wremnants.production.histmaker_tools import (
    aggregate_groups,
    make_quantile_helper,
    scale_to_data,
    write_analysis_output,
)
from wums import logging

parser.add_argument(
    "--csVarsHist", action="store_true", help="Add CS variables to dilepton hist"
)
parser.add_argument("--axes", type=str, nargs="*", default=["mll", "ptll"], help="")
parser.add_argument(
    "--finePtBinning", action="store_true", help="Use fine binning for ptll"
)
parser.add_argument(
    "--useTheoryAgnosticBinning",
    action="store_true",
    help="Use theory agnostic binning (coarser) to produce the results",
)
parser.add_argument(
    "--useDileptonTriggerSelection",
    action="store_true",
    help="Use dilepton trigger selection (default uses the Wlike one, with one triggering muon and odd/even event selection to define its charge, staying agnostic to the other)",
)
parser.add_argument(
    "--insituEffMCFile",
    type=str,
    default=None,
    help="MC efficiency file (.pkl.lz4 from scripts/corrections/make_insitu_effMC.py) "
    "enabling the in-situ muon efficiency systematic histograms (per-category "
    "Chebyshev coefficient variations). Without it, only the effMCprobe_* "
    "histograms needed to produce that file are written.",
)
parser.add_argument(
    "--noAuxiliaryHistograms",
    action="store_true",
    help="Remove auxiliary histograms to save memory (removed by default with --unfolding or --theoryAgnostic)",
)
parser.add_argument(
    "--muonIsolation",
    type=int,
    nargs=2,
    default=[1, 1],
    choices=[-1, 0, 1],
    help="Apply isolation cut to triggering and not-triggering muon (in this order): -1/1 for failing/passing isolation, 0 for skipping it. If using --useDileptonTriggerSelection, then the sorting is based on the muon charge as -/+",
)
parser.add_argument(
    "--flipEventNumberSplitting",
    action="store_true",
    help="Flip even with odd event numbers to consider the positive or negative muon as the W-like muon",
)
parser.add_argument(
    "--useTnpMuonVarForSF",
    action="store_true",
    help="To read efficiency scale factors, use the same muon variables as used to measure them with tag-and-probe (by default the final corrected ones are used)",
)
parser.add_argument(
    "--makeCSQuantileHists",
    action="store_true",
    help="Make hists with fine binned CS variables for producing quantiles",
)
parser.add_argument(
    "--quarkMassCorr",
    nargs="*",
    type=str,
    default=["MiNNLO_Zbb"],
    choices=theory_utils.valid_theory_corrections(),
    help="Apply quark-mass correction generators as additional theory variations.",
)
parser.add_argument(
    "--splitSampleInN",
    type=int,
    default=-1,
    help="Split the sample in N parts, useful for debugging and testing",
)
parser.add_argument(
    "--randomSeedForSplit",
    type=int,
    default=12345,
    help="Random seed for splitting the sample in N parts",
)
parser.add_argument(
    "--jackknifeN",
    type=int,
    default=0,
    help="Number of jackknife samples to use, if > 0, then the sample is split in 2*jackknifeN parts",
)
parser.add_argument(
    "--jackknifeEfficiency",
    type=float,
    default=0.5,
    help="Jackknife efficiency, used to define the size of the sample",
)
parser.add_argument(
    "--randomSeedForJackknife",
    type=int,
    default=12345,
    help="Random seed for jackknifing procedure",
)
parser = parsing.set_parser_default(
    parser, "aggregateGroups", ["Diboson", "Top", "Wtaunu", "Wmunu"]
)
parser = parsing.set_parser_default(parser, "excludeProcs", ["QCD", "DYlowMass"])
parser = parsing.set_parser_default(
    parser, "pt", binning.get_default_ptbins(analysis_label)
)
args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

if args.dxybsVeto > 0 and args.dxybsVeto < args.dxybs:
    raise ValueError("When using together '--dxybsVeto X --dxybs Y' it must be X > Y.")

thisAnalysis = (
    ROOT.wrem.AnalysisType.Dilepton
    if args.useDileptonTriggerSelection
    else ROOT.wrem.AnalysisType.Wlike
)
isoBranch = muon_selections.getIsoBranch(args.isolationDefinition)
era = args.era

datasets = getDatasets(
    maxFiles=args.maxFiles,
    filt=args.filterProcs,
    excl=args.excludeProcs,
    aux=args.auxiliaryProcs,
    nanoVersion="v9",
    base_path=args.dataPath,
    extended="msht20an3lo" not in args.pdfs,
    oneMCfileEveryN=args.oneMCfileEveryN,
    era=era,
)

# dilepton invariant mass cuts
mass_min, mass_max = binning.get_default_mz_window()

ewMassBins = binning.make_bw_binning(mass=91.1535, width=2.4932, initialStep=0.010)

if args.useTheoryAgnosticBinning:
    theoryAgnostic_axes, _ = binning.get_theoryAgnostic_axes(
        ptV_flow=True, absYV_flow=True, wlike=True
    )
    axis_ptV_thag = theoryAgnostic_axes[0]
    dilepton_ptV_binning = axis_ptV_thag.edges
else:
    dilepton_ptV_binning = binning.ptZ_binning if not args.finePtBinning else range(200)

if "yll" in args.axes:
    # use 20 quantiles in case "yll" is used as nominal axis
    edges_yll = binning.yll_20quantiles_binning
    edges_absYll = edges_yll[len(edges_yll) // 2 :]
    axis_yll = hist.axis.Variable(edges_yll, name="yll")
    axis_absYll = hist.axis.Variable(edges_absYll, name="absYll", underflow=False)
else:
    axis_yll = hist.axis.Regular(100, -2.5, 2.5, name="yll")
    axis_absYll = hist.axis.Regular(50, 0.0, 2.5, name="absYll", underflow=False)

# available axes for dilepton validation plots
all_axes = {
    "mll": hist.axis.Variable(
        [
            60,
            70,
            75,
            78,
            80,
            82,
            84,
            85,
            86,
            87,
            88,
            89,
            90,
            91,
            92,
            93,
            94,
            95,
            96,
            97,
            98,
            100,
            102,
            105,
            110,
            120,
        ],
        name="mll",
    ),
    "yll": axis_yll,
    "absYll": axis_absYll,
    "ptll": hist.axis.Variable(dilepton_ptV_binning, name="ptll", underflow=False),
    "etaPlus": hist.axis.Variable([-2.4, -1.2, -0.3, 0.3, 1.2, 2.4], name="etaPlus"),
    "etaMinus": hist.axis.Variable([-2.4, -1.2, -0.3, 0.3, 1.2, 2.4], name="etaMinus"),
    "etaRegionSign": hist.axis.Regular(
        3, 0, 3, name="etaRegionSign", underflow=False, overflow=False
    ),
    "etaRegionRange": hist.axis.Regular(
        3, 0, 3, name="etaRegionRange", underflow=False, overflow=False
    ),
    "absEtaPlus": hist.axis.Regular(8, 0, 2.4, name="absEtaPlus"),
    "absEtaMinus": hist.axis.Regular(8, 0, 2.4, name="absEtaMinus"),
    "etaAbsEta": hist.axis.Variable(
        [
            -2.4,
            -2.0,
            -1.6,
            -1.4,
            -1.2,
            -1.0,
            -0.6,
            0.0,
            0.6,
            1.0,
            1.2,
            1.4,
            1.6,
            2.0,
            2.4,
        ],
        name="etaAbsEta",
    ),
    "etaSum": hist.axis.Regular(12, -4.8, 4.8, name="etaSum"),
    "etaDiff": hist.axis.Variable(
        [-4.8, -1.0, -0.6, -0.2, 0.2, 0.6, 1.0, 4.8], name="etaDiff"
    ),
    "ptPlus": hist.axis.Regular(int(args.pt[0]), args.pt[1], args.pt[2], name="ptPlus"),
    "ptMinus": hist.axis.Regular(
        int(args.pt[0]), args.pt[1], args.pt[2], name="ptMinus"
    ),
    "cosThetaStarll": hist.axis.Regular(
        200 if args.makeCSQuantileHists else 20,
        -1.0,
        1.0,
        name="cosThetaStarll",
        underflow=False,
        overflow=False,
    ),
    "phiStarll": hist.axis.Regular(
        200 if args.makeCSQuantileHists else 20,
        -math.pi,
        math.pi,
        circular=True,
        name="phiStarll",
    ),
    # "charge": hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge") # categorical axes in python bindings always have an overflow bin, so use a regular
    "massVgen": hist.axis.Variable(ewMassBins, name="massVgen"),
    "ewMll": hist.axis.Variable(ewMassBins, name="ewMll"),
    "ewMlly": hist.axis.Variable(ewMassBins, name="ewMlly"),
    "ewLogDeltaM": hist.axis.Regular(100, -10, 4, name="ewLogDeltaM"),
    "firstMuons_abseta0": hist.axis.Regular(
        3, 0.0, 2.4, name="firstMuons_abseta0", underflow=False
    ),
    "secondMuons_eta0": hist.axis.Regular(
        int(args.eta[0]), args.eta[1], args.eta[2], name="secondMuons_eta0"
    ),
    "secondMuons_pt0": hist.axis.Regular(
        int(args.pt[0]), args.pt[1], args.pt[2], name="secondMuons_pt0"
    ),
    "secondMuons_charge0": hist.axis.Regular(
        2, -2.0, 2.0, underflow=False, overflow=False, name="secondMuons_charge0"
    ),
    "ptll_resolution": hist.axis.Regular(1000, -1, 1, name="ptll_resolution"),
}

for a in args.axes:
    if a not in all_axes.keys():
        logger.error(
            f" {a} is not a known axes! Supported axes choices are {list(all_axes.keys())}"
        )

# sepcial axes
axis_corr_eta = hist.axis.Regular(48, -2.4, 2.4, name="corr_eta")
axis_corr_phi = hist.axis.Regular(1, -np.pi, np.pi, circular=True, name="corr_phi")
axis_corr_parms = hist.axis.StrCategory(
    ["A_k", "e_k", "M_k", "M_lambda", "A_phi", "e_phi", "M_phi"],
    name="corr_parms",
)
axis_res_parms = hist.axis.StrCategory(["a", "c", "b"], name="res_parms")

auxiliary_gen_axes = [
    "massVgen",  # preFSR variables
    "ewMll",
    "ewMlly",
    "ewLogDeltaM",  # ew variables
]

# axis order must match cols order (HistoBoost fills positionally): eta, pt, charge
axes_probe_ID = [
    hist.axis.Regular(
        int(args.eta[0]),
        args.eta[1],
        args.eta[2],
        name="eta",
        underflow=False,
        overflow=False,
    ),
    hist.axis.Regular(
        int(args.pt[0]),
        args.pt[1],
        args.pt[2],
        name="pt",
        underflow=False,
        overflow=False,
    ),
    hist.axis.Regular(2, -2.0, 2.0, name="charge", underflow=False, overflow=False),
]
cols_probe_ID = ["probeMuons_eta0", "probeMuons_pt0", "probeMuons_charge0"]

axes_probe_HLT = [
    hist.axis.Regular(
        int(args.eta[0]),
        args.eta[1],
        args.eta[2],
        name="eta",
        underflow=False,
        overflow=False,
    ),
    hist.axis.Regular(
        int(args.pt[0]),
        args.pt[1],
        args.pt[2],
        name="pt",
        underflow=False,
        overflow=False,
    ),
]
cols_probe_HLT = [
    "probeMuons_eta0",
    "probeMuons_pt0",
]

# Probe binning used to measure our own MC efficiency (effMC) for the in-situ
# muon efficiency method: args.eta x args.pt x charge, matching the Chebyshev
# decorrelation grid and probe pt window (see make_insitu_effMC.py).
axes_insitu_effMC = [
    hist.axis.Regular(
        int(args.eta[0]),
        args.eta[1],
        args.eta[2],
        name="eta",
        underflow=False,
        overflow=False,
    ),
    hist.axis.Regular(
        int(args.pt[0]),
        args.pt[1],
        args.pt[2],
        name="pt",
        underflow=False,
        overflow=False,
    ),
    hist.axis.Regular(2, -2.0, 2.0, name="charge", underflow=False, overflow=False),
]
cols_insitu_effMC = ["probeMuons_eta0", "probeMuons_pt0", "probeMuons_charge0"]


nominal_cols = args.axes

if args.csVarsHist:
    # in case CS variables are added to the main histogram, use optimized binning
    # CS variables will be binned in nxn quantiles; quantiles are computed in each bin of args.axes as provided by the quantile_file
    n_quantiles = 8
    all_axes["cosThetaStarll_quantile"] = hist.axis.Regular(
        n_quantiles,
        0,
        1,
        name="cosThetaStarll_quantile",
        underflow=False,
        overflow=False,
    )
    all_axes["phiStarll_quantile"] = hist.axis.Regular(
        n_quantiles,
        0,
        1,
        name="phiStarll_quantile",
        underflow=False,
        overflow=False,
    )

    quantile_file = f"{common.data_dir}/angularCoefficients/mz_dilepton_scetlib_dyturbo_CT18Z_N3p0LL_N2LO_Corr_maxFiles_m1_csQuantiles.hdf5"
    quantile_helper_csVars = make_quantile_helper(
        quantile_file,
        ["cosThetaStarll", "phiStarll"],
        ["ptll", "absYll"],
        name="nominal_csQuantiles",
        processes=["Zmumu_2016PostVFP"],
        n_quantiles=[n_quantiles],
    )

    nominal_cols += ["cosThetaStarll_quantile", "phiStarll_quantile"]

nominal_axes = [all_axes[a] for a in nominal_cols]

if args.unfolding:
    add_helicity_axis = "helicitySig" in args.unfoldingAxes

    if args.unfoldingInclusive:
        cutsmap = {"fiducial": "masswindow"}
    else:
        cutsmap = {
            "pt_min": args.pt[1],
            "pt_max": args.pt[2],
            "abseta_max": args.eta[2],
            "mass_min": mass_min,
            "mass_max": mass_max,
        }

    unfolder_z = unfolding_tools.UnfolderZ(
        reco_axes_edges={a: all_axes[a].edges for a in args.axes},
        unfolding_axes_names=args.unfoldingAxes,
        unfolding_levels=args.unfoldingLevels,
        poi_as_noi=args.poiAsNoi,
        fitresult=args.fitresult,
        cutsmap=cutsmap,
    )

    if not args.poiAsNoi:
        datasets = unfolding_tools.add_out_of_acceptance(datasets, group="Zmumu")

# define helpers
muon_prefiring_helper, muon_prefiring_helper_stat, muon_prefiring_helper_syst = (
    muon_prefiring.make_muon_prefiring_helpers(era=era)
)

if args.skipByHelicityCorrection:
    helicity_smoothing_helpers_procs = {}
else:
    procs = [
        p
        for p, grp in (("W", samples.wprocs), ("Z", samples.zprocs))
        if any(d.name in grp for d in datasets)
    ]
    helicity_smoothing_helpers_procs = (
        theory_corrections.make_helicity_smoothing_helpers(
            args.pdfs, args.theoryCorr, procs=procs
        )
    )

# extra axes which can be used to label tensor_axes
if args.binnedScaleFactors:
    logger.info("Using binned scale factors and uncertainties")
    # might never use it really anymore, but let's warn the user that this is obsolete
    logger.warning(
        "Only SF with no uT dependence are implemented, and the treatment for trigger is like Wlike"
    )
    # add usePseudoSmoothing=True for tests with Asimov
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = (
        muon_efficiencies_binned.make_muon_efficiency_helpers_binned(
            filename=args.sfFile, era=era, max_pt=args.pt[2], is_w_like=True
        )
    )
else:
    logger.info("Using smoothed scale factors and uncertainties")
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = (
        muon_efficiencies_smooth.make_muon_efficiency_helpers_smooth(
            filename=args.sfFile,
            era=era,
            max_pt=args.pt[2],
            what_analysis=thisAnalysis,
            isoEfficiencySmoothing=args.isoEfficiencySmoothing,
            smooth3D=args.smooth3dsf,
            isoDefinition=args.isolationDefinition,
            baseEff_types=["reco", "tracking"],
            trigEff_types=[],
            isoEff_types=[],
        )
    )
logger.info(f"SF file: {args.sfFile}")

# In-situ muon efficiency helper (ID/HLT/Iso floated as Chebyshev polynomials).
# Built only when the MC-efficiency file is provided; the first pass that
# produces that file (via the effMCprobe_* histograms) runs without it.
muon_insitu_efficiency_helper = None
insitu_parameter_labels = None
if args.insituEffMCFile is not None:
    muon_insitu_efficiency_helper, insitu_parameter_labels = (
        muon_efficiencies_insitu.make_muon_insitu_efficiency_helper(
            args.insituEffMCFile, pt_range=(args.pt[1], args.pt[2])
        )
    )
    logger.info(f"In-situ effMC file: {args.insituEffMCFile}")

muon_efficiency_helper_syst_altBkg = {}
for es in common.muonEfficiency_altBkgSyst_effSteps:
    altSFfile = args.sfFile.replace(".root", "_altBkg.root")
    logger.info(f"Additional SF file for alternate syst with {es}: {altSFfile}")
    muon_efficiency_helper_syst_altBkg[es] = (
        muon_efficiencies_smooth.make_muon_efficiency_helpers_smooth_altSyst(
            filename=altSFfile,
            era=era,
            what_analysis=thisAnalysis,
            max_pt=args.pt[2],
            effStep=es,
        )
    )

pileup_helper = pileup.make_pileup_helper(era=era)
vertex_helper = vertex.make_vertex_helper(era=era)

calib_filepaths = common.calib_filepaths
closure_filepaths = common.closure_filepaths
diff_weights_helper = (
    ROOT.wrem.SplinesDifferentialWeightsHelper(calib_filepaths["tflite_file"])
    if (args.muonScaleVariation == "smearingWeightsSplines" or args.validationHists)
    else None
)
(
    mc_jpsi_crctn_helper,
    data_jpsi_crctn_helper,
    mc_jpsi_crctn_unc_helper,
    data_jpsi_crctn_unc_helper,
) = muon_calibration.make_jpsi_crctn_helpers(
    calib_filepaths,
    muon_corr_mc=args.muonCorrMC,
    muon_corr_data=args.muonCorrData,
    scale_var_method=args.muonScaleVariation,
    scale_A=args.scale_A,
    scale_e=args.scale_e,
    scale_M=args.scale_M,
    make_uncertainty_helper=True,
)
z_non_closure_parametrized_helper, z_non_closure_binned_helper = (
    muon_calibration.make_Z_non_closure_helpers(
        args, calib_filepaths, closure_filepaths
    )
)

mc_calibration_helper, data_calibration_helper, calibration_uncertainty_helper = (
    muon_calibration.make_muon_calibration_helpers(args, era=era)
)

closure_unc_helper = muon_calibration.make_closure_uncertainty_helper(
    common.closure_filepaths["parametrized"]
)
closure_unc_helper_A = muon_calibration.make_uniform_closure_uncertainty_helper(
    0, common.correlated_variation_base_size["A"]
)
closure_unc_helper_M = muon_calibration.make_uniform_closure_uncertainty_helper(
    2, common.correlated_variation_base_size["M"]
)

smearing_helper, smearing_uncertainty_helper = (
    (None, None) if args.noSmearing else muon_calibration.make_muon_smearing_helpers()
)

smearinggradhelper = muon_calibration.make_smearing_grad_helper()

bias_helper = muon_calibration.make_muon_bias_helpers(args)

(
    pixel_multiplicity_helper,
    pixel_multiplicity_uncertainty_helper,
    pixel_multiplicity_uncertainty_helper_stat,
) = muon_calibration.make_pixel_multiplicity_helpers(
    reverse_variations=args.reweightPixelMultiplicity
)

if args.nToysMC > 0:
    seed_data = 2 * args.randomSeedForToys
    seed_mc = 2 * args.randomSeedForToys + 1
    toy_helper_data = ROOT.wrem.ToyHelper(
        args.nToysMC, seed_data, 1, ROOT.ROOT.GetThreadPoolSize()
    )
    toy_helper_mc = ROOT.wrem.ToyHelper(
        args.nToysMC,
        seed_mc,
        args.varianceScalingForToys,
        ROOT.ROOT.GetThreadPoolSize(),
    )
    axis_toys = hist.axis.Integer(
        0, args.nToysMC, underflow=False, overflow=False, name="toys"
    )
if args.splitSampleInN > 1:
    seed_mc_split = 2 * args.randomSeedForSplit + 2
    rand_helper_mc = ROOT.wrem.RandomUniformHelper(
        args.splitSampleInN, seed_mc_split, ROOT.ROOT.GetThreadPoolSize()
    )
    axis_split = hist.axis.Integer(
        0,
        args.splitSampleInN,
        underflow=False,
        overflow=False,
        name="sample_split",
    )
if args.jackknifeN > 0:
    seed_mc_jackknife = 2 * args.randomSeedForJackknife + 1
    jackknife_helper = ROOT.wrem.JackknifeHelper(
        args.jackknifeN,
        args.jackknifeEfficiency,
        seed_mc_jackknife,
        ROOT.ROOT.GetThreadPoolSize(),
    )
    axis_jackknife = hist.axis.Integer(
        0, args.jackknifeN, underflow=False, overflow=False, name="jackknife_sample"
    )

procs_v = [d.name for d in datasets if d.name in samples.vprocs]
theory_corrs = [*args.theoryCorr, *args.ewTheoryCorr]
corr_helpers = theory_corrections.load_corr_helpers(procs_v, theory_corrs)
if args.quarkMassCorr:
    procs_z = [d.name for d in datasets if d.name in samples.zprocs]
    corr_helpers_quark_mass = theory_corrections.load_corr_helpers(
        procs_z, args.quarkMassCorr
    )
    for proc, helper_map in corr_helpers_quark_mass.items():
        corr_helpers.setdefault(proc, {})
        corr_helpers[proc].update(helper_map)


def build_graph(df, dataset):
    logger.info(f"build graph for dataset: {dataset.name}")
    results = []
    isW = dataset.name in samples.wprocs
    isZ = dataset.name in samples.zprocs
    isWorZ = isW or isZ

    if isWorZ and dataset.name[0] in helicity_smoothing_helpers_procs.keys():
        helicity_smoothing_helpers = helicity_smoothing_helpers_procs[dataset.name[0]]
    else:
        helicity_smoothing_helpers = {}

    cvh_helper = data_calibration_helper if dataset.is_data else mc_calibration_helper
    jpsi_helper = data_jpsi_crctn_helper if dataset.is_data else mc_jpsi_crctn_helper
    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    df = df.DefinePerSample("unity", "1.0")
    df = df.Define(
        "isEvenEvent", f"event % 2 {'!=' if args.flipEventNumberSplitting else '=='} 0"
    )

    if args.nToysMC > 0:
        if dataset.is_data:
            df = df.Define("toyIdxs", toy_helper_data, ["rdfslot_"])
        else:
            df = df.Define("toyIdxs", toy_helper_mc, ["rdfslot_"])

    if args.splitSampleInN > 1 and not dataset.is_data:
        df = df.Define("sample_n", rand_helper_mc, ["rdfslot_"])

    if args.jackknifeN > 0 and not dataset.is_data:
        df = df.Define("jackknife_sample", jackknife_helper, ["rdfslot_"])

    weightsum = df.SumAndCount("weight")

    axes = nominal_axes
    cols = nominal_cols

    if args.addRunAxis and dataset.is_data:
        run_edges = binning.run_edges
        axes = [
            *axes,
            hist.axis.Variable(
                run_edges + 0.5, name="run", underflow=False, overflow=False
            ),
        ]
        cols = [*cols, "run"]

    if args.unfolding and dataset.group == "Zmumu":
        print(f"name = {dataset.name}; group = {dataset.group}")
        df = unfolder_z.add_gen_histograms(
            args,
            df,
            results,
            dataset,
            corr_helpers,
            helicity_smoothing_helpers=helicity_smoothing_helpers,
        )

        if not unfolder_z.poi_as_noi:
            axes = [
                *nominal_axes,
                *unfolder_z.unfolding_axes[unfolder_z.unfolding_levels[-1]],
            ]
            cols = [
                *nominal_cols,
                *unfolder_z.unfolding_cols[unfolder_z.unfolding_levels[-1]],
            ]

    if args.xnormOnly:
        return results, weightsum

    if not args.noAuxiliaryHistograms and isZ and len(auxiliary_gen_axes):
        # gen level variables before selection
        df_gen = df
        df_gen = df_gen.DefinePerSample("exp_weight", "1.0")
        df_gen = theory_corrections.define_theory_weights_and_corrs(
            df_gen,
            dataset.name,
            corr_helpers,
            args,
            helicity_smoothing_helpers=helicity_smoothing_helpers,
        )

        for obs in auxiliary_gen_axes:
            results.append(
                df_gen.HistoBoost(
                    f"gen_{obs}", [all_axes[obs]], [obs, "nominal_weight"]
                )
            )
            systematics.add_theory_hists(
                results,
                df_gen,
                args,
                dataset.name,
                corr_helpers,
                helicity_smoothing_helpers,
                [all_axes[obs]],
                [obs],
                base_name=f"gen_{obs}",
                for_wmass=False,
            )

    df = df.Filter(muon_selections.hlt_string(era))

    df = muon_selections.veto_electrons(df)
    df = muon_selections.apply_met_filters(df)

    df = muon_calibration.define_corrected_muons(
        df, cvh_helper, jpsi_helper, args, dataset, smearing_helper, bias_helper
    )

    df = muon_selections.select_veto_muons(
        df,
        nMuons=2,
        etaCut=args.vetoRecoEta,
        staPtCut=args.vetoRecoStaPt,
        dxybsCut=args.dxybsVeto if args.dxybsVeto > 0 else args.dxybs,
    )
    isoThreshold = args.isolationThreshold
    passIsoBoth = args.muonIsolation[0] + args.muonIsolation[1] == 2
    df = muon_selections.select_good_muons(
        df,
        args.pt[1],
        args.pt[2],
        dataset.group,
        nMuons=1,
        use_trackerMuons=args.trackerMuons,
        use_isolation=passIsoBoth,
        isoBranch=isoBranch,
        isoThreshold=isoThreshold,
        requirePixelHits=args.requirePixelHits,
        dxybsCut=args.dxybs,
        condition=">=",
    )

    df = muon_selections.define_two_muons(df, dilepton=args.useDileptonTriggerSelection)

    # iso cut applied here, if requested, because it needs the definition of firstMuons and secondMuons from muon_selections.define_trigger_muons
    if not passIsoBoth and (args.muonIsolation[0] or args.muonIsolation[1]):
        df = muon_selections.apply_iso_muons(
            df,
            args.muonIsolation[0],
            args.muonIsolation[1],
            isoBranch,
            isoThreshold,
            name_first="firstMuons",
            name_second="secondMuons",
        )

    df = muon_selections.select_z_candidate(
        df, mass_min, mass_max, name_first="firstMuons", name_second="secondMuons"
    )

    df = muon_selections.select_standalone_muons(
        df, dataset, args.trackerMuons, "firstMuons"
    )
    df = muon_selections.select_standalone_muons(
        df, dataset, args.trackerMuons, "secondMuons"
    )

    if args.useDileptonTriggerSelection:
        df = muon_selections.apply_triggermatching_muon(
            df, dataset, "firstMuons", "secondMuons", era=era
        )
        df = df.Alias("muonsMinus_pt0", "firstMuons_pt0")
        df = df.Alias("muonsPlus_pt0", "secondMuons_pt0")
        df = df.Alias("muonsMinus_eta0", "firstMuons_eta0")
        df = df.Alias("muonsPlus_eta0", "secondMuons_eta0")
        df = df.Alias("muonsMinus_mom4", "firstMuons_mom4")
        df = df.Alias("muonsPlus_mom4", "secondMuons_mom4")
    else:
        df = muon_selections.apply_triggermatching_muon(
            df, dataset, "firstMuons", era=era
        )
        df = df.Define("firstMuon_isNegative", "firstMuons_charge0 == -1")
        df = df.Define(
            "muonsMinus_pt0", "firstMuon_isNegative ? firstMuons_pt0 : secondMuons_pt0"
        )
        df = df.Define(
            "muonsPlus_pt0", "firstMuon_isNegative ? secondMuons_pt0 : firstMuons_pt0"
        )
        df = df.Define(
            "muonsMinus_eta0",
            "firstMuon_isNegative ? firstMuons_eta0 : secondMuons_eta0",
        )
        df = df.Define(
            "muonsPlus_eta0",
            "firstMuon_isNegative ? secondMuons_eta0 : firstMuons_eta0",
        )
        df = df.Define(
            "muonsMinus_mom4",
            "firstMuon_isNegative ? firstMuons_mom4 : secondMuons_mom4",
        )
        df = df.Define(
            "muonsPlus_mom4",
            "firstMuon_isNegative ? secondMuons_mom4 : firstMuons_mom4",
        )

    useTnpMuonVarForSF = args.useTnpMuonVarForSF
    # in principle these are only needed for MC,
    # but one may want to compare tnp and corrected variables also for data
    if useTnpMuonVarForSF:
        df = df.Define("firstMuons_tnpPt0", "Muon_pt[firstMuons][0]")
        df = df.Define("firstMuons_tnpEta0", "Muon_eta[firstMuons][0]")
        df = df.Define("firstMuons_tnpCharge0", "Muon_charge[firstMuons][0]")
        df = df.Define("secondMuons_tnpPt0", "Muon_pt[secondMuons][0]")
        df = df.Define("secondMuons_tnpEta0", "Muon_eta[secondMuons][0]")
        df = df.Define("secondMuons_tnpCharge0", "Muon_charge[secondMuons][0]")
    else:
        df = df.Alias("firstMuons_tnpPt0", "firstMuons_pt0")
        df = df.Alias("firstMuons_tnpEta0", "firstMuons_eta0")
        df = df.Alias("firstMuons_tnpCharge0", "firstMuons_charge0")
        df = df.Alias("secondMuons_tnpPt0", "secondMuons_pt0")
        df = df.Alias("secondMuons_tnpEta0", "secondMuons_eta0")
        df = df.Alias("secondMuons_tnpCharge0", "secondMuons_charge0")

    df = df.Define("ptll", "ll_mom4.pt()")
    df = df.Define("yll", "ll_mom4.Rapidity()")
    df = df.Define("absYll", "std::fabs(yll)")
    # "renaming" to write out corresponding axis
    df = df.Alias("ptMinus", "muonsMinus_pt0")
    df = df.Alias("ptPlus", "muonsPlus_pt0")
    df = df.Alias("etaMinus", "muonsMinus_eta0")
    df = df.Alias("etaPlus", "muonsPlus_eta0")
    df = df.Define("absEtaMinus", "std::fabs(etaMinus)")
    df = df.Define("absEtaPlus", "std::fabs(etaPlus)")
    df = df.Define("etaAbsEta", "absEtaMinus > absEtaPlus ? etaMinus : etaPlus")

    df = df.Define(
        "etaRegionRange",
        "(std::abs(muonsPlus_eta0) > 0.9) + (std::abs(muonsMinus_eta0) > 0.9)",
    )  # eta region: 0: barrel-barrel, 1: endcap-barrel, 2: endcap-endcap
    df = df.Define(
        "etaRegionSign", "(muonsPlus_eta0 > 0) + (muonsMinus_eta0 > 0)"
    )  # eta region: 0: both muons in negative eta, 1: one muon in negative eta, 2: both muons in positive eta

    df = df.Define("etaSum", "muonsPlus_eta0 + muonsMinus_eta0")
    df = df.Define("etaDiff", "muonsPlus_eta0 - muonsMinus_eta0")  # plus - minus

    df = df.Define(
        "csSineCosThetaPhill",
        "wrem::csSineCosThetaPhi(muonsPlus_mom4, muonsMinus_mom4)",
    )
    df = df.Define("cosThetaStarll", "csSineCosThetaPhill.costheta")
    df = df.Define("phiStarll", "csSineCosThetaPhill.phi()")

    if args.csVarsHist:
        for c, h, a in (
            (
                "phiStarll_quantile",
                quantile_helper_csVars[0],
                ["phiStarll", "ptll", "absYll"],
            ),
            (
                "cosThetaStarll_quantile",
                quantile_helper_csVars[1],
                ["cosThetaStarll", "phiStarll", "ptll", "absYll"],
            ),
        ):
            if [a for a in h.axes.name] != a:
                raise RuntimeError(
                    f"Invalid helper axes: {[a for a in h.axes.name]} != {a}"
                )

            df = df.Define(c, h, a)

    # TODO might need to add an explicit cut on firstMuons_pt0 in case nominal pt range
    # extends below 26 GeV e.g. for calibration test purposes
    df = df.Define("firstMuons_abseta0", "std::fabs(firstMuons_eta0)")

    if dataset.is_data:
        df = df.DefinePerSample("nominal_weight", "1.0")
        cvhName = "cvh"
    else:
        cvhName = "cvhideal"

    axis_eta = hist.axis.Regular(int(args.eta[0]), args.eta[1], args.eta[2], name="eta")
    axis_pt = hist.axis.Regular(int(args.pt[0]), args.pt[1], args.pt[2], name="pt")
    axis_charge = binning.axis_charge
    axis_nvalidpixel = hist.axis.Integer(0, 10, name="nvalidpixel")

    df = df.Define(
        f"firstMuons_{cvhName}NValidPixelHits0",
        f"Muon_{cvhName}NValidPixelHits[firstMuons][0]",
    )
    df = df.Define(
        f"secondMuons_{cvhName}NValidPixelHits0",
        f"Muon_{cvhName}NValidPixelHits[secondMuons][0]",
    )

    logger.debug(f"Define weights and store nominal histograms")

    df = df.Define("goodFirstMuons", "goodMuons && firstMuons")
    df = df.Define("goodSecondMuons", "goodMuons && secondMuons")

    df = df.Define("firstMuons_passID0", "Sum(goodFirstMuons) == 1")
    df = df.Define("secondMuons_passID0", "Sum(goodSecondMuons) == 1")

    if not dataset.is_data:
        df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
        df = df.Define("weight_vtx", vertex_helper, ["GenVtx_z", "Pileup_nTrueInt"])

        if era == "2016PostVFP":
            df = df.Define(
                "weight_newMuonPrefiringSF",
                muon_prefiring_helper,
                [
                    "Muon_correctedEta",
                    "Muon_correctedPt",
                    "Muon_correctedPhi",
                    "Muon_correctedCharge",
                    "Muon_looseId",
                ],
            )
            weight_expr = (
                "weight_pu*weight_newMuonPrefiringSF*L1PreFiringWeight_ECAL_Nom"
            )
        else:
            weight_expr = (
                "weight_pu*L1PreFiringWeight_Muon_Nom*L1PreFiringWeight_ECAL_Nom"
            )

        if not args.noVertexWeight:
            weight_expr += "*weight_vtx"

        muonVarsForSF = [
            "tnpPt0",
            "tnpEta0",
            "SApt0",
            "SAeta0",
            # "tnpUT0",
            "tnpCharge0",
            # "passIso0",
        ]
        # if args.useDileptonTriggerSelection:
        # muonVarsForSF.append("passTrigger0")
        # careful, first all firstMuon variables, then all secondMuon
        columnsForSF = [
            f"{t}Muons_{v}" for t in ["first", "second"] for v in muonVarsForSF
        ]

        df = muon_selections.define_muon_uT_variable(
            df,
            isWorZ,
            smooth3dsf=args.smooth3dsf,
            colNamePrefix="firstMuons",
            addWithTnpMuonVar=useTnpMuonVarForSF,
        )
        df = muon_selections.define_muon_uT_variable(
            df,
            isWorZ,
            smooth3dsf=args.smooth3dsf,
            colNamePrefix="secondMuons",
            addWithTnpMuonVar=useTnpMuonVarForSF,
        )
        # # ut is defined in muon_selections.define_muon_uT_variable
        # if not useTnpMuonVarForSF:
        #     df = df.Alias("firstMuons_tnpUT0", "firstMuons_uT0")
        #     df = df.Alias("secondMuons_tnpUT0", "secondMuons_uT0")

        # if not args.smooth3dsf:
        #     columnsForSF.remove("firstMuons_tnpUT0")
        #     columnsForSF.remove("secondMuons_tnpUT0")

        # reco + tracking scale factors stay from the external measurement
        # (ID/HLT/Iso are determined in-situ instead). Applied to the nominal
        # weight for all categories.
        if not args.noScaleFactors:
            df = df.Define(
                "weight_fullMuonSF_withTrackingReco",
                muon_efficiency_helper,
                columnsForSF,
            )
            weight_expr += "*weight_fullMuonSF_withTrackingReco"

        # prepare inputs for pixel multiplicity helpers.
        # NOTE: the pixel-multiplicity correction has a (nonTriggering, triggering)
        # category axis (the shared wrem::TriggerCat enum). In this dilepton
        # selection both muons can trigger, so the category is assigned
        # POSITIONALLY (first muon -> triggering slot, second muon ->
        # nonTriggering slot); it does NOT reflect the actual per-muon trigger
        # decision. The RVec element order is irrelevant (the helpers loop over
        # muons), only the per-muon (kinematics <-> category) pairing matters.
        df = df.DefinePerSample(
            "MuonFirstSecond_triggerCat",
            "ROOT::VecOps::RVec<wrem::TriggerCat>{wrem::TriggerCat::triggering, wrem::TriggerCat::nonTriggering}",
        )
        df = df.Define(
            "MuonFirstSecond_eta",
            "ROOT::VecOps::RVec<float>{firstMuons_eta0, secondMuons_eta0}",
        )
        df = df.Define(
            "MuonFirstSecond_pt",
            "ROOT::VecOps::RVec<float>{firstMuons_pt0, secondMuons_pt0}",
        )
        df = df.Define(
            "MuonFirstSecond_charge",
            "ROOT::VecOps::RVec<int>{firstMuons_charge0, secondMuons_charge0}",
        )
        df = df.Define(
            f"MuonFirstSecond_{cvhName}NValidPixelHits",
            f"ROOT::VecOps::RVec<int>{{firstMuons_{cvhName}NValidPixelHits0, secondMuons_{cvhName}NValidPixelHits0}}",
        )

        pixel_multiplicity_cols = [
            "MuonFirstSecond_triggerCat",
            "MuonFirstSecond_eta",
            "MuonFirstSecond_pt",
            "MuonFirstSecond_charge",
            f"MuonFirstSecond_{cvhName}NValidPixelHits",
        ]

        if args.reweightPixelMultiplicity:
            df = df.Define(
                "weight_pixel_multiplicity",
                pixel_multiplicity_helper,
                pixel_multiplicity_cols,
            )
            weight_expr += "*weight_pixel_multiplicity"

        logger.debug(f"Experimental weight defined: {weight_expr}")
        df = df.Define("exp_weight", weight_expr)
        df = theory_corrections.define_theory_weights_and_corrs(
            df,
            dataset.name,
            corr_helpers,
            args,
            helicity_smoothing_helpers=helicity_smoothing_helpers,
        )

        df = df.Define(
            "Mupluscor_mom4",
            "firstMuons_charge0 == 1 ? firstMuons_mom4 : secondMuons_mom4",
        )
        df = df.Define(
            "Muminuscor_mom4",
            "firstMuons_charge0 == -1 ? firstMuons_mom4 : secondMuons_mom4",
        )

        df = df.Define(
            "parmgrads_k7",
            "wrem::parmgrads_k7_t(Mupluscor_mom4, Muminuscor_mom4, nominal_weight)",
        )
        df = df.Define(
            "parmgradsres_k3",
            smearinggradhelper,
            ["Mupluscor_mom4", "Muminuscor_mom4", "nominal_weight"],
        )

        df = df.Define(
            "etav", "std::array<double, 2>{Mupluscor_mom4.eta(), Muminuscor_mom4.eta()}"
        )
        df = df.Define(
            "phiv", "std::array<double, 2>{Mupluscor_mom4.phi(), Muminuscor_mom4.phi()}"
        )
        df = df.Define(
            "ptv", "std::array<double, 2>{Mupluscor_mom4.pt(), Muminuscor_mom4.pt()}"
        )

    df = df.Define("firstMuons_tag0", "firstMuons_passID0 && firstMuons_passTrigger0")
    df = df.Define(
        "secondMuons_tag0", "secondMuons_passID0 && secondMuons_passTrigger0"
    )

    # Select events with exactly one muon that passes ID and HLT
    df_1HLT = df.Filter("firstMuons_tag0 != secondMuons_tag0")

    df_1HLT = df_1HLT.Define("probeMuons", "firstMuons_tag0 ? secondMuons : firstMuons")
    df_1HLT = df_1HLT.Define(
        "probeMuons_passID0",
        "firstMuons_tag0 ? secondMuons_passID0 : firstMuons_passID0",
    )
    df_1HLT = df_1HLT.Define(
        "probeMuons_pt0", "firstMuons_tag0 ? secondMuons_pt0 : firstMuons_pt0"
    )
    df_1HLT = df_1HLT.Define(
        "probeMuons_eta0", "firstMuons_tag0 ? secondMuons_eta0 : firstMuons_eta0"
    )
    df_1HLT = df_1HLT.Define(
        "probeMuons_charge0",
        "firstMuons_tag0 ? secondMuons_charge0 : firstMuons_charge0",
    )
    # tag is the other leg (always passes ID & HLT)
    df_1HLT = df_1HLT.Define(
        "tagMuons_pt0", "firstMuons_tag0 ? firstMuons_pt0 : secondMuons_pt0"
    )
    df_1HLT = df_1HLT.Define(
        "tagMuons_eta0", "firstMuons_tag0 ? firstMuons_eta0 : secondMuons_eta0"
    )
    df_1HLT = df_1HLT.Define(
        "tagMuons_charge0", "firstMuons_tag0 ? firstMuons_charge0 : secondMuons_charge0"
    )

    df_1HLT_passID = df_1HLT.Filter("probeMuons_passID0 == 1")
    df_1HLT_failID = df_1HLT.Filter("probeMuons_passID0 == 0")

    df_2HLT = df.Filter("firstMuons_tag0 && secondMuons_tag0")
    # the two muons are indistinguishable and valid tag muons, take randomly one or the other as probe
    df_2HLT = df_2HLT.Define("probeMuons", "isEvenEvent ? firstMuons : secondMuons")
    df_2HLT = df_2HLT.Define(
        "probeMuons_passIso0", f"{isoBranch}[probeMuons][0] < {isoThreshold}"
    )
    df_2HLT = df_2HLT.Define(
        "probeMuons_pt0", "isEvenEvent ? firstMuons_pt0 : secondMuons_pt0"
    )
    df_2HLT = df_2HLT.Define(
        "probeMuons_eta0", "isEvenEvent ? firstMuons_eta0 : secondMuons_eta0"
    )
    df_2HLT = df_2HLT.Define(
        "probeMuons_charge0", "isEvenEvent ? firstMuons_charge0 : secondMuons_charge0"
    )
    # tag is the other leg (also passes ID & HLT in the 2HLT category)
    df_2HLT = df_2HLT.Define(
        "tagMuons_pt0", "isEvenEvent ? secondMuons_pt0 : firstMuons_pt0"
    )
    df_2HLT = df_2HLT.Define(
        "tagMuons_eta0", "isEvenEvent ? secondMuons_eta0 : firstMuons_eta0"
    )
    df_2HLT = df_2HLT.Define(
        "tagMuons_charge0", "isEvenEvent ? secondMuons_charge0 : firstMuons_charge0"
    )

    df_2HLT_passIso = df_2HLT.Filter("probeMuons_passIso0 == 1")
    df_2HLT_failIso = df_2HLT.Filter("probeMuons_passIso0 == 0")

    dfs = {
        "nominal": {"df": df_2HLT_passIso, "axes": axes, "cols": cols},
        "failIso": {
            "df": df_2HLT_failIso,
            "axes": axes_probe_HLT,
            "cols": cols_probe_HLT,
        },
        "failHLT": {
            "df": df_1HLT_passID,
            "axes": axes_probe_HLT,
            "cols": cols_probe_HLT,
        },
        "failID": {"df": df_1HLT_failID, "axes": axes_probe_ID, "cols": cols_probe_ID},
    }

    # In-situ effMC measurement (MC only): emit the probe (eta, pt, charge)
    # spectrum of each of the 4 categories. make_insitu_effMC.py recombines
    # these into per-step MC efficiencies (idip/trigger/iso), e.g.
    #   eff_iso     = nominal / (nominal + failIso)
    #   eff_trigger = (nominal + failIso) / (nominal + failIso + failHLT)
    #   eff_idip    = (nominal + failIso + failHLT)
    #                 / (nominal + failIso + failHLT + failID)
    if not dataset.is_data:
        for _cat, _info in dfs.items():
            results.append(
                _info["df"].HistoBoost(
                    f"effMCprobe_{_cat}",
                    axes_insitu_effMC,
                    [*cols_insitu_effMC, "nominal_weight"],
                )
            )

    for channel, info in dfs.items():
        df = info["df"]
        axes = info["axes"]
        cols = info["cols"]

        if dataset.is_data:
            if args.nToysMC > 0:
                axes = [*axes, axis_toys]
                cols = [*cols, "toyIdxs"]

            results.append(df.HistoBoost(channel, axes, cols))
        else:
            results.append(
                df.HistoBoost(
                    "weight",
                    [hist.axis.Regular(100, -2, 2)],
                    ["nominal_weight"],
                    storage=hist.storage.Double(),
                )
            )

            if args.nToysMC > 0 or args.splitSampleInN > 1 or args.jackknifeN > 1:
                results.append(
                    df.HistoBoost(f"{channel}_asimov", axes, [*cols, "nominal_weight"])
                )
            if args.nToysMC > 0:
                axes = [*axes, axis_toys]
                cols = [*cols, "toyIdxs"]
            if args.splitSampleInN > 1:
                axes = [*axes, axis_split]
                cols = [*cols, "sample_n"]
            if args.jackknifeN > 1:
                axes = [*axes, axis_jackknife]
                cols = [*cols, "jackknife_sample"]

            results.append(df.HistoBoost(channel, axes, [*cols, "nominal_weight"]))

            if isZ:
                # theory agnostic stuff
                theoryAgnostic_axes, theoryAgnostic_cols = (
                    binning.get_theoryAgnostic_axes(
                        ptV_bins=[],
                        absYV_bins=[],
                        ptV_flow=True,
                        absYV_flow=True,
                        wlike=True,
                    )
                )
                axis_helicity = binning.axis_helicity_multidim

                df_theory_agnostic = theoryAgnostic_tools.define_helicity_weights(
                    df, is_z=True
                )
                noiAsPoiHistName = common.hist_name(
                    channel, syst="yieldsTheoryAgnostic"
                )
                logger.debug(
                    f"Creating special histogram '{noiAsPoiHistName}' for theory agnostic to treat POIs as NOIs"
                )
                results.append(
                    df_theory_agnostic.HistoBoost(
                        noiAsPoiHistName,
                        [*axes, *theoryAgnostic_axes],
                        [*cols, *theoryAgnostic_cols, "nominal_weight_helicity"],
                        tensor_axes=[axis_helicity],
                    )
                )

            if args.unfolding and args.poiAsNoi and dataset.group == "Zmumu":
                unfolder_z.add_poi_as_noi_histograms(
                    df, results, axes, cols, channel=channel
                )

        # histograms for corrections/uncertainties for pixel hit multiplicity

        # valid-pixel-hit spectra per muon leg (inputs to make_pixel_corrections.py).
        # first/second are the negative/positive muon; the pixel correction's
        # triggering/nonTriggering slots are filled positionally from first/second
        # (see make_pixel_corrections.py and the note above).
        hNValidPixelHitsFirst = df.HistoBoost(
            f"{channel}_hNValidPixelHitsFirst",
            [axis_eta, axis_pt, axis_charge, axis_nvalidpixel],
            [
                "firstMuons_eta0",
                "firstMuons_pt0",
                "firstMuons_charge0",
                f"firstMuons_{cvhName}NValidPixelHits0",
                "nominal_weight",
            ],
        )
        results.append(hNValidPixelHitsFirst)

        hNValidPixelHitsSecond = df.HistoBoost(
            f"{channel}_hNValidPixelHitsSecond",
            [axis_eta, axis_pt, axis_charge, axis_nvalidpixel],
            [
                "secondMuons_eta0",
                "secondMuons_pt0",
                "secondMuons_charge0",
                f"secondMuons_{cvhName}NValidPixelHits0",
                "nominal_weight",
            ],
        )
        results.append(hNValidPixelHitsSecond)

        if args.makeCSQuantileHists:
            results.append(
                df.HistoBoost(
                    f"{channel}_csQuantiles",
                    [
                        all_axes[o]
                        for o in ["ptll", "absYll", "phiStarll", "cosThetaStarll"]
                    ],
                    ["ptll", "absYll", "phiStarll", "cosThetaStarll"],
                )
            )

        if not args.noAuxiliaryHistograms:
            for obs in [
                ["ptll", "yll"],
                "mll",
                "cosThetaStarll",
                "phiStarll",
                "etaPlus",
                "etaMinus",
                "ptPlus",
                "ptMinus",
            ]:
                if isinstance(obs, str):
                    obs = [obs]
                obs_name = "_".join([channel, *obs])
                obs_axes = [all_axes[o] for o in obs]

                if dataset.is_data:
                    results.append(df.HistoBoost(obs_name, obs_axes, obs))
                else:
                    results.append(
                        df.HistoBoost(obs_name, obs_axes, [*obs, "nominal_weight"])
                    )
                    if isWorZ and not args.onlyMainHistograms:
                        df = systematics.add_theory_hists(
                            results,
                            df,
                            args,
                            dataset.name,
                            corr_helpers,
                            helicity_smoothing_helpers,
                            obs_axes,
                            obs,
                            base_name=obs_name,
                            for_wmass=False,
                        )

        if not args.noAuxiliaryHistograms and isZ:
            # gen level variables
            for obs in auxiliary_gen_axes:
                results.append(
                    df.HistoBoost(
                        f"{channel}_{obs}", [all_axes[obs]], [obs, "nominal_weight"]
                    )
                )
                if not args.onlyMainHistograms:
                    df = systematics.add_theory_hists(
                        results,
                        df,
                        args,
                        dataset.name,
                        corr_helpers,
                        helicity_smoothing_helpers,
                        [all_axes[obs]],
                        [obs],
                        base_name=f"{channel}_{obs}",
                        for_wmass=False,
                    )

        if not dataset.is_data:
            parmgrad_axes = [*axes[:-1], axis_corr_eta, axis_corr_phi]
            parmgrad_cols = [*cols[:-1], "etav", "phiv", "parmgrads_k7"]

            hparmgrads = df.HistoBoost(
                f"{channel}_hparmgrads",
                parmgrad_axes,
                parmgrad_cols,
                tensor_axes=[axis_corr_parms],
            )
            results.append(hparmgrads)

            parmgradres_axes = parmgrad_axes
            parmgradres_cols = [*cols[:-1], "etav", "phiv", "parmgradsres_k3"]

            hparmgradsres = df.HistoBoost(
                f"{channel}_hparmgradsres",
                parmgradres_axes,
                parmgradres_cols,
                tensor_axes=[axis_res_parms],
            )
            results.append(hparmgradsres)

        if not dataset.is_data and not args.onlyMainHistograms:

            # reco/tracking efficiency uncertainties stay from the external SF
            # measurement (ID/HLT/Iso are handled in-situ instead)
            df = systematics.add_muon_efficiency_unc_hists(
                results,
                df,
                muon_efficiency_helper_stat,
                muon_efficiency_helper_syst,
                axes,
                cols,
                what_analysis=thisAnalysis,
                smooth3D=args.smooth3dsf,
                base_name=channel,
            )

            if muon_insitu_efficiency_helper is not None:
                # in-situ ID/HLT/Iso efficiency: per-category Chebyshev
                # coefficient variations (unconstrained nuisances in the fit)
                df = systematics.add_muon_insitu_efficiency_hists(
                    results,
                    df,
                    muon_insitu_efficiency_helper,
                    axes,
                    cols,
                    category=channel,
                    base_name=channel,
                )

            df = systematics.add_L1Prefire_unc_hists(
                results,
                df,
                axes,
                cols,
                helper_stat=muon_prefiring_helper_stat,
                helper_syst=muon_prefiring_helper_syst,
                base_name=channel,
            )

            if isWorZ:

                df = systematics.add_theory_hists(
                    results,
                    df,
                    args,
                    dataset.name,
                    corr_helpers,
                    helicity_smoothing_helpers,
                    axes,
                    cols,
                    for_wmass=False,
                    base_name=channel,
                )

                reco_sel = "vetoMuonsPre"
                require_prompt = "tau" not in dataset.name
                df = muon_calibration.define_genFiltered_recoMuonSel(
                    df, reco_sel, require_prompt
                )
                reco_sel_GF = muon_calibration.getColName_genFiltered_recoMuonSel(
                    reco_sel, require_prompt
                )
                df = muon_calibration.define_matched_gen_muons_kinematics(
                    df, reco_sel_GF
                )
                df = muon_calibration.calculate_matched_gen_muon_kinematics(
                    df, reco_sel_GF
                )
                df = muon_calibration.define_matched_reco_muon_kinematics(
                    df, reco_sel_GF
                )

                ####################################################
                # nuisances from the muon momemtum scale calibration
                if args.muonCorrData in ["massfit", "lbl_massfit"]:
                    input_kinematics = [
                        f"{reco_sel_GF}_recoPt",
                        f"{reco_sel_GF}_recoEta",
                        f"{reco_sel_GF}_recoCharge",
                        f"{reco_sel_GF}_genPt",
                        f"{reco_sel_GF}_genEta",
                        f"{reco_sel_GF}_genCharge",
                    ]
                    if diff_weights_helper:
                        df = df.Define(
                            f"{reco_sel_GF}_response_weight",
                            diff_weights_helper,
                            [*input_kinematics],
                        )
                        input_kinematics.append(f"{reco_sel_GF}_response_weight")

                    # muon scale variation from stats. uncertainty on the jpsi massfit
                    df = df.Define(
                        "nominal_muonScaleSyst_responseWeights_tensor",
                        data_jpsi_crctn_unc_helper,
                        [*input_kinematics, "nominal_weight"],
                    )
                    muonScaleSyst_responseWeights = df.HistoBoost(
                        f"{channel}_muonScaleSyst_responseWeights",
                        axes,
                        [*cols, "nominal_muonScaleSyst_responseWeights_tensor"],
                        tensor_axes=data_jpsi_crctn_unc_helper.tensor_axes,
                        storage=hist.storage.Double(),
                    )
                    results.append(muonScaleSyst_responseWeights)

                    df = muon_calibration.add_resolution_uncertainty(
                        df,
                        axes,
                        results,
                        cols,
                        smearing_uncertainty_helper,
                        reco_sel_GF,
                        base_name=channel,
                    )

                    # add pixel multiplicity uncertainties
                    df = df.Define(
                        "nominal_pixelMultiplicitySyst_tensor",
                        pixel_multiplicity_uncertainty_helper,
                        [*pixel_multiplicity_cols, "nominal_weight"],
                    )
                    hist_pixelMultiplicitySyst = df.HistoBoost(
                        f"{channel}_pixelMultiplicitySyst",
                        axes,
                        [*cols, "nominal_pixelMultiplicitySyst_tensor"],
                        tensor_axes=pixel_multiplicity_uncertainty_helper.tensor_axes,
                        storage=hist.storage.Double(),
                    )
                    results.append(hist_pixelMultiplicitySyst)

                    if args.pixelMultiplicityStat:
                        df = df.Define(
                            "nominal_pixelMultiplicityStat_tensor",
                            pixel_multiplicity_uncertainty_helper_stat,
                            [*pixel_multiplicity_cols, "nominal_weight"],
                        )
                        hist_pixelMultiplicityStat = df.HistoBoost(
                            f"{channel}_pixelMultiplicityStat",
                            axes,
                            [*cols, "nominal_pixelMultiplicityStat_tensor"],
                            tensor_axes=pixel_multiplicity_uncertainty_helper_stat.tensor_axes,
                            storage=hist.storage.Double(),
                        )
                        results.append(hist_pixelMultiplicityStat)

                    if args.nonClosureScheme in ["A-M-separated", "A-only"]:
                        # add the ad-hoc Z non-closure nuisances from the jpsi massfit to muon scale unc
                        df = df.DefinePerSample("AFlag", "0x01")
                        df = df.Define(
                            "Z_non_closure_parametrized_A",
                            z_non_closure_parametrized_helper,
                            [*input_kinematics, "nominal_weight", "AFlag"],
                        )
                        hist_Z_non_closure_parametrized_A = df.HistoBoost(
                            f"{channel}_Z_non_closure_parametrized_A",
                            axes,
                            [*cols, "Z_non_closure_parametrized_A"],
                            tensor_axes=z_non_closure_parametrized_helper.tensor_axes,
                            storage=hist.storage.Double(),
                        )
                        results.append(hist_Z_non_closure_parametrized_A)

                    if args.nonClosureScheme in [
                        "A-M-separated",
                        "binned-plus-M",
                        "M-only",
                    ]:
                        df = df.DefinePerSample("MFlag", "0x04")
                        df = df.Define(
                            "Z_non_closure_parametrized_M",
                            z_non_closure_parametrized_helper,
                            [*input_kinematics, "nominal_weight", "MFlag"],
                        )
                        hist_Z_non_closure_parametrized_M = df.HistoBoost(
                            f"{channel}_Z_non_closure_parametrized_M",
                            axes,
                            [*cols, "Z_non_closure_parametrized_M"],
                            tensor_axes=z_non_closure_parametrized_helper.tensor_axes,
                            storage=hist.storage.Double(),
                        )
                        results.append(hist_Z_non_closure_parametrized_M)

                    if args.nonClosureScheme == "A-M-combined":
                        df = df.DefinePerSample("AMFlag", "0x01 | 0x04")
                        df = df.Define(
                            "Z_non_closure_parametrized",
                            z_non_closure_parametrized_helper,
                            [*input_kinematics, "nominal_weight", "AMFlag"],
                        )
                        hist_Z_non_closure_parametrized = df.HistoBoost(
                            (
                                f"{channel}_Z_non_closure_parametrized_gaus"
                                if args.muonScaleVariation == "smearingWeightsGaus"
                                else f"{channel}_Z_non_closure_parametrized"
                            ),
                            axes,
                            [*cols, "Z_non_closure_parametrized"],
                            tensor_axes=z_non_closure_parametrized_helper.tensor_axes,
                            storage=hist.storage.Double(),
                        )
                        results.append(hist_Z_non_closure_parametrized)

                    # extra uncertainties from non-closure stats
                    df = df.Define(
                        "muonScaleClosSyst_responseWeights_tensor_splines",
                        closure_unc_helper,
                        [*input_kinematics, "nominal_weight"],
                    )
                    nominal_muonScaleClosSyst_responseWeights = df.HistoBoost(
                        f"{channel}_muonScaleClosSyst_responseWeights",
                        axes,
                        [*cols, "muonScaleClosSyst_responseWeights_tensor_splines"],
                        tensor_axes=closure_unc_helper.tensor_axes,
                        storage=hist.storage.Double(),
                    )
                    results.append(nominal_muonScaleClosSyst_responseWeights)

                    # extra uncertainties for A (fully correlated)
                    df = df.Define(
                        "muonScaleClosASyst_responseWeights_tensor_splines",
                        closure_unc_helper_A,
                        [*input_kinematics, "nominal_weight"],
                    )
                    nominal_muonScaleClosASyst_responseWeights = df.HistoBoost(
                        f"{channel}_muonScaleClosASyst_responseWeights",
                        axes,
                        [*cols, "muonScaleClosASyst_responseWeights_tensor_splines"],
                        tensor_axes=closure_unc_helper_A.tensor_axes,
                        storage=hist.storage.Double(),
                    )
                    results.append(nominal_muonScaleClosASyst_responseWeights)

                    # extra uncertainties for M (fully correlated)
                    df = df.Define(
                        "muonScaleClosMSyst_responseWeights_tensor_splines",
                        closure_unc_helper_M,
                        [*input_kinematics, "nominal_weight"],
                    )
                    nominal_muonScaleClosMSyst_responseWeights = df.HistoBoost(
                        f"{channel}_muonScaleClosMSyst_responseWeights",
                        axes,
                        [*cols, "muonScaleClosMSyst_responseWeights_tensor_splines"],
                        tensor_axes=closure_unc_helper_M.tensor_axes,
                        storage=hist.storage.Double(),
                    )
                    results.append(nominal_muonScaleClosMSyst_responseWeights)

                ####################################################

                # Don't think it makes sense to apply the mass weights to scale leptons from tau decays
                if not "tau" in dataset.name:
                    systematics.add_muonscale_hist(
                        results,
                        df,
                        args.muonCorrEtaBins,
                        args.muonCorrMag,
                        isW,
                        axes,
                        cols,
                        muon_eta="firstMuons_eta0",
                        base_name=channel,
                    )  ## FIXME: what muon to choose ?

    if hasattr(dataset, "out_of_acceptance"):
        # Rename dataset to not overwrite the original one
        dataset.name = dataset.name + "OOA"

    return results, weightsum


logger.debug(f"Datasets are {[d.name for d in datasets]}")
resultdict = narf.build_and_run(datasets[::-1], build_graph)

if not args.noScaleToData:
    scale_to_data(resultdict)
    aggregate_groups(datasets, resultdict, args.aggregateGroups)

write_analysis_output(
    resultdict, f"{os.path.basename(__file__).replace('py', 'hdf5')}", args
)
