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
    "--insituEffMCFile",
    type=str,
    default=None,
    help="MC efficiency file (.pkl.lz4 from scripts/corrections/make_insitu_effMC.py) "
    "enabling the in-situ muon efficiency systematic histograms (per-category "
    "Chebyshev coefficient variations). Without it, only the effMCprobe_* "
    "histograms needed to produce that file are written.",
)
parser.add_argument(
    "--insituSFFile",
    type=str,
    default=None,
    help="Accumulated central in-situ Chebyshev coefficients (theta_central) "
    ".pkl.lz4 from scripts/corrections/make_insitu_effSF.py, for the iterative "
    "fit. The MC nominal is reweighted by the post-fit central SF before all "
    "category templates and syst tensors are filled, and the in-situ variations "
    "are re-linearised around it. Default None = iteration 0 (theta_central=0, "
    "== SF=1, current behaviour).",
)
parser.add_argument(
    "--makeInsituEffMC",
    action="store_true",
    help="Emit the effMCprobe_* probe spectra (both tag-and-probe legs per "
    "2HLT event) used by scripts/corrections/make_insitu_effMC.py to compute "
    "the in-situ MC efficiencies. Only needed for that one-time pre-step.",
)
parser.add_argument(
    "--makeUTQuantileHists",
    action="store_true",
    help="Emit fine-binned recoUT input histograms in (eta, pt, [charge]) for "
    "failIso/failHLT, used to compute per-cell uT quantile edges in a follow-up "
    "make_quantile_helper pass.",
)
parser.add_argument(
    "--utQuantileFile",
    type=str,
    default=None,
    help="Path to the histmaker output produced with --makeUTQuantileHists; when "
    "set, the failIso/failHLT recoUT fit axis is replaced by a per-(eta, pt, "
    "[charge]) 7-quantile axis.",
)
parser.add_argument(
    "--noAuxiliaryHistograms",
    action="store_true",
    help="Remove auxiliary histograms to save memory (removed by default with --unfolding or --theoryAgnostic)",
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

thisAnalysis = ROOT.wrem.AnalysisType.Dilepton

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

axis_eta = hist.axis.Regular(
    int(args.eta[0]),
    args.eta[1],
    args.eta[2],
    name="eta",
    underflow=False,
    overflow=False,
)
axis_pt = hist.axis.Regular(
    int(args.pt[0]), args.pt[1], args.pt[2], name="pt", underflow=False, overflow=False
)
axis_charge = binning.axis_charge

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

# failID has falling probe-pt stats above ~45 GeV (few/empty bins on the 1 GeV grid
# at high pt × forward eta). Use a variable-width pt axis: 1 GeV up to 38 GeV, then
# coarser bins toward the configured upper pt cut args.pt[2] (last edge follows the
# cut, so --pt 44 26 70 widens the top bin to [50, 70]; args.pt[2] == 60 reproduces
# the original axis). This same axis is reused for the in-situ effMC grid
# (axes_insitu_effMC) so effMC and the fit share binning; the Chebyshev x̃ window is
# set separately from (args.pt[1], args.pt[2]).
axis_pt_failID = hist.axis.Variable(
    [
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        38,
        40,
        42,
        45,
        50,
        int(round(args.pt[2])),
    ],
    name="pt",
    underflow=False,
    overflow=False,
)

# reco-uT axis for the failHLT and failIso fit histograms. Range matches the
# SMP-23-002 smoothSF3D_uTm30to100 convention. Coarsened in the tails relative
# to the original 10-bin grid because the Z-pT distribution makes |uT| > 20 GeV
# rare -> avoid empty bins in the high-|uT| corners.
axis_recoUT = hist.axis.Variable(
    [-30, -10, -5, 0, 5, 10, 20, 40, 100],
    name="recoUT",
    underflow=False,
    overflow=False,
)
# Fine-binned recoUT axis used only as the *target axis* of the quantile
# pre-step (--makeUTQuantileHists); a Regular 1 GeV grid gives enough
# resolution for the per-cell CDF interpolation that make_quantile_helper does.
axis_recoUT_fine = hist.axis.Regular(
    130,
    -30,
    100,
    name="recoUT",
    underflow=False,
    overflow=False,
)
# Quantile axes used as the *fit axis* when --utQuantileFile is set: events are
# mapped to a quantile fraction in [0, 1) by the per-cell CDF helper, then binned
# here. failHLT (trigger) uses 10 quantiles for finer reco-uT resolution (the
# trigger uT/pt coefficients are the under-constrained ones); failIso stays at 7.
n_ut_quantiles_failHLT = 10
n_ut_quantiles_failIso = 7
axis_recoUT_quantile_failHLT = hist.axis.Regular(
    n_ut_quantiles_failHLT,
    0,
    1,
    name="recoUT_quantile",
    underflow=False,
    overflow=False,
)
axis_recoUT_quantile_failIso = hist.axis.Regular(
    n_ut_quantiles_failIso,
    0,
    1,
    name="recoUT_quantile",
    underflow=False,
    overflow=False,
)
# Gen-uT axis for the effMC probe spectra: same coarse 8-bin grid as recoUT.
axis_genUT = hist.axis.Variable(
    [-30, -10, -5, 0, 5, 10, 20, 40, 100],
    name="genUT",
    underflow=False,
    overflow=False,
)

# axis order must match cols order (HistoBoost fills positionally)
# failID: eta, pt(variable), charge -> IDIP has no uT dependence
axes_probe_ID = [axis_eta, axis_pt_failID, axis_charge]
cols_probe_ID = ["probeMuons_eta0", "probeMuons_pt0", "probeMuons_charge0"]

if args.utQuantileFile is not None:
    # Per-(eta, pt, [charge]) uT quantiles: each cell has its own 7 equal-
    # population uT bins. The recoUT_quantile column is defined per-category
    # in the dfs block via make_quantile_helper.
    axes_probe_failHLT = [
        axis_eta,
        axis_pt_failID,
        axis_recoUT_quantile_failHLT,
        axis_charge,
    ]
    cols_probe_failHLT = [
        "probeMuons_eta0",
        "probeMuons_pt0",
        "probeMuons_recoUT_quantile",
        "probeMuons_charge0",
    ]
    axes_probe_failIso = [axis_eta, axis_pt_failID, axis_recoUT_quantile_failIso]
    cols_probe_failIso = [
        "probeMuons_eta0",
        "probeMuons_pt0",
        "probeMuons_recoUT_quantile",
    ]
else:
    # failHLT: eta, pt(variable), recoUT, charge -> Trigger is 2D (pt, uT), charge-dep.
    # Variable pt axis (same as failID) avoids empty bins at high probe pt.
    axes_probe_failHLT = [axis_eta, axis_pt_failID, axis_recoUT, axis_charge]
    cols_probe_failHLT = [
        "probeMuons_eta0",
        "probeMuons_pt0",
        "probeMuons_recoUT0",
        "probeMuons_charge0",
    ]
    # failIso: eta, pt(variable), recoUT -> Iso is 2D (pt, uT), charge-inclusive.
    axes_probe_failIso = [axis_eta, axis_pt_failID, axis_recoUT]
    cols_probe_failIso = [
        "probeMuons_eta0",
        "probeMuons_pt0",
        "probeMuons_recoUT0",
    ]

# Probe binning used to measure our own MC efficiency (effMC) for the in-situ
# muon efficiency method. The genUT axis is needed for HLT and Iso effMC
# (sliced/summed appropriately in make_insitu_effMC.py).
axes_insitu_effMC = [
    hist.axis.Regular(
        int(args.eta[0]),
        args.eta[1],
        args.eta[2],
        name="eta",
        underflow=False,
        overflow=False,
    ),
    # reuse the fail-category variable pt axis so effMC and the fit share binning
    # (the helper looks up effMC via boost-hist axis.index(), variable-safe)
    axis_pt_failID,
    hist.axis.Regular(2, -2.0, 2.0, name="charge", underflow=False, overflow=False),
    axis_genUT,
]
cols_insitu_effMC = [
    "probeMuons_eta0",
    "probeMuons_pt0",
    "probeMuons_charge0",
    "probeMuons_tnpUT0",
]

# Per-cell uT quantile helpers (made once from the pre-step file). Each helper
# returns a quantile fraction in [0, 1) for the event's recoUT, given the
# event's (eta, pt[variable], [charge]) cell.
quantile_helper_failIso = None
quantile_helper_failHLT = None
if args.utQuantileFile is not None:
    quantile_helper_failIso = make_quantile_helper(
        args.utQuantileFile,
        ["recoUT"],
        ["eta", "pt"],
        name="utQuantileInput_failIso",
        processes=["Zmumu_2016PostVFP"],
    )[0]
    quantile_helper_failHLT = make_quantile_helper(
        args.utQuantileFile,
        ["recoUT"],
        ["eta", "pt", "charge"],
        name="utQuantileInput_failHLT",
        processes=["Zmumu_2016PostVFP"],
    )[0]


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
muon_insitu_central_helper = None
insitu_parameter_labels = None
# The in-situ central reweight folds W(theta_central) into nominal_weight, which
# also feeds the effMCprobe_* spectra. effMC must stay iteration-independent, so
# the effMC pre-step and an iteration reweight cannot run together.
if args.makeInsituEffMC and args.insituSFFile is not None:
    raise ValueError(
        "--makeInsituEffMC and --insituSFFile are mutually exclusive: the central "
        "reweight would distort the (iteration-independent) effMC. Compute effMC at "
        "iteration 0 (no --insituSFFile), then reuse it for all iterations."
    )
if args.insituEffMCFile is not None:
    (
        muon_insitu_efficiency_helper,
        muon_insitu_central_helper,
        insitu_parameter_labels,
    ) = muon_efficiencies_insitu.make_muon_insitu_efficiency_helper(
        args.insituEffMCFile,
        pt_range=(args.pt[1], args.pt[2]),
        central_sf_file=args.insituSFFile,
    )
    logger.info(f"In-situ effMC file: {args.insituEffMCFile}")
    if args.insituSFFile is not None:
        logger.info(f"In-situ central SF file (iteration): {args.insituSFFile}")

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
    df = muon_selections.select_good_muons(
        df,
        args.pt[1],
        args.pt[2],
        nMuons=1,
        use_trackerMuons=args.trackerMuons,
        use_isolation=False,
        isoBranch=isoBranch,
        isoThreshold=isoThreshold,
        requirePixelHits=args.requirePixelHits,
        dxybsCut=args.dxybs,
        condition=">=",
    )

    df = muon_selections.define_two_muons(df, dilepton=True)

    df = muon_selections.select_z_candidate(
        df, mass_min, mass_max, name_first="firstMuons", name_second="secondMuons"
    )

    df = muon_selections.select_standalone_muons(
        df, dataset, args.trackerMuons, "firstMuons"
    )
    df = muon_selections.select_standalone_muons(
        df, dataset, args.trackerMuons, "secondMuons"
    )

    df = muon_selections.apply_triggermatching_muon(
        df, dataset, "firstMuons", "secondMuons", era=era
    )
    df = df.Alias("muonsMinus_pt0", "firstMuons_pt0")
    df = df.Alias("muonsPlus_pt0", "secondMuons_pt0")
    df = df.Alias("muonsMinus_eta0", "firstMuons_eta0")
    df = df.Alias("muonsPlus_eta0", "secondMuons_eta0")
    df = df.Alias("muonsMinus_mom4", "firstMuons_mom4")
    df = df.Alias("muonsPlus_mom4", "secondMuons_mom4")

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
    # reco Z phi for the reco-uT fit axis on failHLT/failIso histograms.
    df = df.Define("phill", "ll_mom4.Phi()")
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
            "tnpCharge0",
        ]
        # careful, first all firstMuon variables, then all secondMuon
        columnsForSF = [
            f"{t}Muons_{v}" for t in ["first", "second"] for v in muonVarsForSF
        ]

        # Always define uT (gen-level projection of the boson onto the muon
        # pT direction). Needed by the in-situ helper for the (pT, uT) HLT/Iso
        # Chebyshev terms (matches SMP-23-002 convention). Force smooth3dsf=True
        # so the real uT is computed and not the dummy 0.0f fallback.
        df = muon_selections.define_muon_uT_variable(
            df,
            isWorZ,
            smooth3dsf=True,
            colNamePrefix="firstMuons",
            addWithTnpMuonVar=useTnpMuonVarForSF,
        )
        df = muon_selections.define_muon_uT_variable(
            df,
            isWorZ,
            smooth3dsf=True,
            colNamePrefix="secondMuons",
            addWithTnpMuonVar=useTnpMuonVarForSF,
        )
        if not useTnpMuonVarForSF:
            # When useTnpMuonVarForSF is False, only _uT0 is defined by
            # define_muon_uT_variable. Alias so _tnpUT0 always exists.
            df = df.Alias("firstMuons_tnpUT0", "firstMuons_uT0")
            df = df.Alias("secondMuons_tnpUT0", "secondMuons_uT0")

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

        # Iterative in-situ SF (iteration >= 1): fold the post-fit central
        # reweight W(theta_central) into the *experimental* weight here -- exactly
        # like the reco/tracking SF above -- so it is part of nominal_weight
        # BEFORE define_theory_weights_and_corrs builds the pdf/scetlib/EW weight
        # tensors. Those tensors bake in nominal_weight, so they must see W too;
        # otherwise syst/nominal (= logk) is distorted by 1/W per bin. Compute the
        # per-event probe/tag + 4-category assignment on the parent df (same logic
        # as the tag-and-probe split below) with insituW_-prefixed columns to
        # avoid colliding with the split's probeMuons_*/firstMuons_tag0. Events
        # with no tag muon (dropped from every category downstream) get W = 1.
        # Gated on --insituSFFile (iteration 0 has theta_central=0 -> W=1, so we
        # skip and stay bit-identical to the un-iterated run).
        if args.insituSFFile is not None and muon_insitu_central_helper is not None:
            _ci = systematics.insitu_category_index
            df = df.Define(
                "insituW_firstTag", "firstMuons_passID0 && firstMuons_passTrigger0"
            )
            df = df.Define(
                "insituW_secondTag", "secondMuons_passID0 && secondMuons_passTrigger0"
            )
            df = df.Define(
                "insituW_firstIso", f"{isoBranch}[firstMuons][0] < {isoThreshold}"
            )
            df = df.Define(
                "insituW_secondIso", f"{isoBranch}[secondMuons][0] < {isoThreshold}"
            )
            # probe leg: 2HLT -> random (isEvenEvent); 1HLT -> the non-tag leg
            df = df.Define(
                "insituW_probeIsFirst",
                "(insituW_firstTag && insituW_secondTag) ? isEvenEvent "
                ": static_cast<bool>(!insituW_firstTag)",
            )
            for _v in ["pt0", "eta0", "charge0", "tnpUT0"]:
                df = df.Define(
                    f"insituW_probe_{_v}",
                    f"insituW_probeIsFirst ? firstMuons_{_v} : secondMuons_{_v}",
                )
                df = df.Define(
                    f"insituW_tag_{_v}",
                    f"insituW_probeIsFirst ? secondMuons_{_v} : firstMuons_{_v}",
                )
            # 4-category index (matches systematics.insitu_category_index)
            df = df.Define(
                "insituW_cat",
                "static_cast<int>((insituW_firstTag && insituW_secondTag) ? "
                "((insituW_probeIsFirst ? insituW_firstIso : insituW_secondIso) ? "
                f"{_ci['nominal']} : {_ci['failIso']}) : "
                "((insituW_probeIsFirst ? firstMuons_passID0 : secondMuons_passID0) ? "
                f"{_ci['failHLT']} : {_ci['failID']}))",
            )
            df = df.Define(
                "insituCentralW_raw",
                muon_insitu_central_helper,
                [
                    "insituW_probe_pt0",
                    "insituW_probe_eta0",
                    "insituW_probe_charge0",
                    "insituW_probe_tnpUT0",
                    "insituW_tag_pt0",
                    "insituW_tag_eta0",
                    "insituW_tag_charge0",
                    "insituW_tag_tnpUT0",
                    "insituW_cat",
                ],
            )
            # events with no tag muon are dropped from every category -> W = 1
            df = df.Define(
                "insituCentralW",
                "(insituW_firstTag + insituW_secondTag) >= 1 ? insituCentralW_raw : 1.0",
            )
            weight_expr += "*insituCentralW"

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
    df_1HLT = df_1HLT.Define(
        "probeMuons_phi0",
        "firstMuons_tag0 ? secondMuons_phi0 : firstMuons_phi0",
    )
    # reco-uT for the probe leg: project probe onto the reco Z (mu1+mu2)
    # direction. Same physical quantity as the gen-level uT used in
    # SMP-23-002 (since p_T_had ~= -p_T(Z)), but computable from reco
    # variables only -> available on Data as well as MC.
    df_1HLT = df_1HLT.Define(
        "probeMuons_recoUT0",
        "static_cast<double>(wrem::zqtproj0_boson(probeMuons_pt0, probeMuons_phi0, ptll, phill))",
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
    if not dataset.is_data:
        # gen-level uT used by the in-situ helper (matches SMP-23-002 SF
        # parameterisation; only MC events are reweighted)
        df_1HLT = df_1HLT.Define(
            "probeMuons_tnpUT0",
            "firstMuons_tag0 ? secondMuons_tnpUT0 : firstMuons_tnpUT0",
        )
        df_1HLT = df_1HLT.Define(
            "tagMuons_tnpUT0",
            "firstMuons_tag0 ? firstMuons_tnpUT0 : secondMuons_tnpUT0",
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
    df_2HLT = df_2HLT.Define(
        "probeMuons_phi0", "isEvenEvent ? firstMuons_phi0 : secondMuons_phi0"
    )
    df_2HLT = df_2HLT.Define(
        "probeMuons_recoUT0",
        "static_cast<double>(wrem::zqtproj0_boson(probeMuons_pt0, probeMuons_phi0, ptll, phill))",
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
    if not dataset.is_data:
        df_2HLT = df_2HLT.Define(
            "probeMuons_tnpUT0",
            "isEvenEvent ? firstMuons_tnpUT0 : secondMuons_tnpUT0",
        )
        df_2HLT = df_2HLT.Define(
            "tagMuons_tnpUT0",
            "isEvenEvent ? secondMuons_tnpUT0 : firstMuons_tnpUT0",
        )

    df_2HLT_passIso = df_2HLT.Filter("probeMuons_passIso0 == 1")
    df_2HLT_failIso = df_2HLT.Filter("probeMuons_passIso0 == 0")

    # Per-cell uT quantile column: helper trained on Zmumu MC (eta, pt, [q])
    # CDF of recoUT; applied to *every* event (Data + MC) so that the failIso
    # / failHLT fit histograms are coherent across processes. Same pattern as
    # the cs-quantile helper (see cosThetaStarll_quantile usage upstream),
    # which works directly because its inputs (ptll, phiStarll, ...) are
    # intrinsically `double`. Our probeMuons_*0 columns are `float`/`int`, so
    # wrap them in `_q` columns cast to `double` to match the narf helper
    # signature (one-line, side effect free — the original columns keep their
    # types and downstream consumers are unaffected).
    # The quantile helper (narf make_hist_helper) does *unchecked* axis
    # indexing, and its recoUT axis is flow-less (Regular(130,-30,100)). recoUT
    # is unbounded (high-pt Z tails exceed +/-the window), so clamp it strictly
    # inside [-30, 100) before the lookup to avoid an out-of-bounds segfault.
    # eta/pt are already in range by the good-muon selection.
    _recoUT_clamp = (
        "std::min(std::max(probeMuons_recoUT0, "
        f"{axis_recoUT_fine.edges[0]}), {axis_recoUT_fine.edges[-1]} - 1e-4)"
    )
    if quantile_helper_failIso is not None:
        df_2HLT_failIso = (
            df_2HLT_failIso.Define("probeMuons_recoUT0_q", _recoUT_clamp)
            .Define("probeMuons_eta0_q", "static_cast<double>(probeMuons_eta0)")
            .Define("probeMuons_pt0_q", "static_cast<double>(probeMuons_pt0)")
            .Define(
                "probeMuons_recoUT_quantile",
                quantile_helper_failIso,
                [
                    "probeMuons_recoUT0_q",
                    "probeMuons_eta0_q",
                    "probeMuons_pt0_q",
                ],
            )
        )
    if quantile_helper_failHLT is not None:
        df_1HLT_passID = (
            df_1HLT_passID.Define("probeMuons_recoUT0_q", _recoUT_clamp)
            .Define("probeMuons_eta0_q", "static_cast<double>(probeMuons_eta0)")
            .Define("probeMuons_pt0_q", "static_cast<double>(probeMuons_pt0)")
            .Define("probeMuons_charge0_q", "static_cast<double>(probeMuons_charge0)")
            .Define(
                "probeMuons_recoUT_quantile",
                quantile_helper_failHLT,
                [
                    "probeMuons_recoUT0_q",
                    "probeMuons_eta0_q",
                    "probeMuons_pt0_q",
                    "probeMuons_charge0_q",
                ],
            )
        )

    dfs = {
        "nominal": {"df": df_2HLT_passIso, "axes": axes, "cols": cols},
        "failIso": {
            "df": df_2HLT_failIso,
            "axes": axes_probe_failIso,
            "cols": cols_probe_failIso,
        },
        "failHLT": {
            "df": df_1HLT_passID,
            "axes": axes_probe_failHLT,
            "cols": cols_probe_failHLT,
        },
        "failID": {"df": df_1HLT_failID, "axes": axes_probe_ID, "cols": cols_probe_ID},
    }

    if args.makeUTQuantileHists and not dataset.is_data:
        # Fine-binned recoUT input hists used by make_quantile_helper to derive
        # per-cell (eta, pt[variable], [charge]) uT CDFs. Order of axes must be
        # (target, *dependent_axes) -- this is what make_quantile_helper expects.
        results.append(
            df_2HLT_failIso.HistoBoost(
                "utQuantileInput_failIso",
                [axis_recoUT_fine, axis_eta, axis_pt_failID],
                [
                    "probeMuons_recoUT0",
                    "probeMuons_eta0",
                    "probeMuons_pt0",
                    "nominal_weight",
                ],
            )
        )
        results.append(
            df_1HLT_passID.HistoBoost(
                "utQuantileInput_failHLT",
                [axis_recoUT_fine, axis_eta, axis_pt_failID, axis_charge],
                [
                    "probeMuons_recoUT0",
                    "probeMuons_eta0",
                    "probeMuons_pt0",
                    "probeMuons_charge0",
                    "nominal_weight",
                ],
            )
        )

    # In-situ effMC measurement (MC only): emit the probe (eta, pt, charge, uT)
    # spectrum of each of the 4 categories. make_insitu_effMC.py recombines
    # these into per-step MC efficiencies (idip/trigger/iso), e.g.
    #   eff_iso     = nominal / (nominal + failIso)
    #   eff_trigger = (nominal + failIso) / (nominal + failIso + failHLT)
    #   eff_idip    = (nominal + failIso + failHLT)
    #                 / (nominal + failIso + failHLT + failID)
    #
    # Dilepton tag-and-probe counting: a 2HLT event provides TWO valid
    # (tag, probe) pairs -- either leg can be the probe (both pass ID &
    # trigger) -- whereas a 1HLT event provides a single valid probe (the
    # non-tag leg; the triggering leg is not a valid probe since its tag would
    # have to pass the trigger, which the other leg fails). We therefore fill
    # BOTH legs of each 2HLT event into the nominal/failIso histograms, each
    # routed by ITS OWN isolation at ITS OWN (eta, pt, charge, uT); 1HLT events
    # fill the single probe leg. The two 2HLT legs always have opposite charge,
    # so they never share a (eta, pt, charge, uT) bin -> clean per-bin stats,
    # and no double-counting factor is needed in the producer.
    if args.makeInsituEffMC and not dataset.is_data:
        # both legs of the 2HLT events, split into iso-pass (nominal) and
        # iso-fail (failIso) probe sets via per-leg isolation masks.
        df_2HLT_eff = df_2HLT.Define(
            "firstMuons_passIso0", f"{isoBranch}[firstMuons][0] < {isoThreshold}"
        )
        df_2HLT_eff = df_2HLT_eff.Define(
            "secondMuons_passIso0", f"{isoBranch}[secondMuons][0] < {isoThreshold}"
        )
        df_2HLT_eff = df_2HLT_eff.Define(
            "effProbe2HLT_isoMask",
            "ROOT::VecOps::RVec<bool>{firstMuons_passIso0, secondMuons_passIso0}",
        )
        df_2HLT_eff = df_2HLT_eff.Define(
            "effProbe2HLT_failIsoMask",
            "ROOT::VecOps::RVec<bool>{!firstMuons_passIso0, !secondMuons_passIso0}",
        )
        df_2HLT_eff = df_2HLT_eff.Define(
            "effProbe2HLT_eta",
            "ROOT::VecOps::RVec<float>{firstMuons_eta0, secondMuons_eta0}",
        )
        df_2HLT_eff = df_2HLT_eff.Define(
            "effProbe2HLT_pt",
            "ROOT::VecOps::RVec<float>{firstMuons_pt0, secondMuons_pt0}",
        )
        df_2HLT_eff = df_2HLT_eff.Define(
            "effProbe2HLT_charge",
            "ROOT::VecOps::RVec<int>{firstMuons_charge0, secondMuons_charge0}",
        )
        df_2HLT_eff = df_2HLT_eff.Define(
            "effProbe2HLT_uT",
            "ROOT::VecOps::RVec<float>{firstMuons_tnpUT0, secondMuons_tnpUT0}",
        )
        for _cat, _mask in (
            ("nominal", "effProbe2HLT_isoMask"),
            ("failIso", "effProbe2HLT_failIsoMask"),
        ):
            for _v in ("eta", "pt", "charge", "uT"):
                df_2HLT_eff = df_2HLT_eff.Define(
                    f"effProbe_{_cat}_{_v}", f"effProbe2HLT_{_v}[{_mask}]"
                )
            results.append(
                df_2HLT_eff.HistoBoost(
                    f"effMCprobe_{_cat}",
                    axes_insitu_effMC,
                    [
                        f"effProbe_{_cat}_eta",
                        f"effProbe_{_cat}_pt",
                        f"effProbe_{_cat}_charge",
                        f"effProbe_{_cat}_uT",
                        "nominal_weight",
                    ],
                )
            )
        # 1HLT categories: single valid probe (the non-tag leg), scalar fill.
        for _cat, _df in (("failHLT", df_1HLT_passID), ("failID", df_1HLT_failID)):
            results.append(
                _df.HistoBoost(
                    f"effMCprobe_{_cat}",
                    axes_insitu_effMC,
                    [*cols_insitu_effMC, "nominal_weight"],
                )
            )

    # NOTE: the in-situ central reweight W(theta_central) is applied upstream as
    # part of nominal_weight (see the insituCentralW block before exp_weight), so
    # by here every category's templates AND theory/syst tensors already carry it.
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

            for es in common.muonEfficiency_altBkgSyst_effSteps:
                df = systematics.add_muon_efficiency_unc_hists_altBkg(
                    results,
                    df,
                    muon_efficiency_helper_syst_altBkg[es],
                    axes,
                    cols,
                    what_analysis=thisAnalysis,
                    step=es,
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
