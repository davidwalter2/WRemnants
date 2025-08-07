import copy
import os
import pprint

import h5py
import hist
import numpy as np

import rabbit
import rabbit.io_tools
from rabbit.tensorwriter import TensorWriter
from utilities import common, parsing
from utilities.io_tools import input_tools
from wremnants import theory_corrections
from wremnants.datasets.datagroups import Datagroups
from wums import boostHistHelpers as hh
from wums import logging


class AlphaSTheoryFitTW(TensorWriter):
    """
    Tensor writer for the alphaS fit direct to theory predictions.
    Will format histograms correctly before adding them as data, processes, variations.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.logger = logging.child_logger(__name__)
        self.ref = {}

    def set_reference(self, channel, h, lumi=16800, scale=1.0, postOp=None):
        self.ref[channel] = {
            "h": h,
            "lumi": lumi,
            "scale": scale,
            "postOp": postOp,
            "ptV_name": self.get_ptV_axis_name(h),
            "absYV_name": self.get_absYV_axis_name(h),
            "chargeV_name": self.get_charge_axis_name(h),
            "ptV_bins": h.axes[self.get_ptV_axis_name(h)].edges,
            "absYV_bins": h.axes[self.get_absYV_axis_name(h)].edges,
        }
        self.logger.debug(f"Initialized channel {channel} with parameters")
        self.logger.debug(pprint.pformat(self.ref[channel]))

    def add_systematic(
        self,
        h,
        name,
        process,
        channel,
        rebin_pt=True,
        rebin_y=True,
        normalize=True,
        apply_postOp=True,
        format=True,
        **kwargs,
    ):
        """
        Systematic variations on a process. Formatting is applied.
        """

        if format:
            if isinstance(h, (list, tuple)):
                h[0] = self.format(
                    h[0],
                    channel,
                    process,
                    rebin_pt=rebin_pt,
                    rebin_y=rebin_y,
                    normalize=normalize,
                    apply_postOp=apply_postOp,
                )
                h[1] = self.format(
                    h[1],
                    channel,
                    process,
                    rebin_pt=rebin_pt,
                    rebin_y=rebin_y,
                    normalize=normalize,
                    apply_postOp=apply_postOp,
                )

            elif kwargs.get("mirror"):
                h = self.format(
                    h,
                    channel,
                    process,
                    rebin_pt=rebin_pt,
                    rebin_y=rebin_y,
                    normalize=normalize,
                    apply_postOp=apply_postOp,
                )

        super().add_systematic(h, name, process, channel, **kwargs)

    def add_scale_systematic(
        self,
        h,
        name,
        process,
        channel,
        rebin_pt=True,
        rebin_y=True,
        normalize=True,
        apply_postOp=True,
        **kwargs,
    ):
        """
        Systematic variations that are computed as
            var_up/var_nominal * nominal
            var_down/var_nominal * nominal (or mirror)
        where nominal is the process in this channel.
        First applies the formatting, then computes the variation.
        """

        if not kwargs.get("mirror"):
            h[0] = self.format(
                h[0],
                channel,
                process,
                rebin_pt=rebin_pt,
                rebin_y=rebin_y,
                normalize=normalize,
                apply_postOp=apply_postOp,
            )
            h[1] = self.format(
                h[1],
                channel,
                process,
                rebin_pt=rebin_pt,
                rebin_y=rebin_y,
                normalize=normalize,
                apply_postOp=apply_postOp,
            )
            h[2] = self.format(
                h[2],
                channel,
                process,
                rebin_pt=rebin_pt,
                rebin_y=rebin_y,
                normalize=normalize,
                apply_postOp=apply_postOp,
            )

            hup = hh.divideHists(h[0], h[2])
            hdown = hh.divideHists(h[1], h[2])
            print(hup)
            print(self.ref[channel][process])
            hup = hh.multiplyHists(hup, self.ref[channel][process])
            hdown = hh.multiplyHists(hdown, self.ref[channel][process])

            super().add_systematic([hup, hdown], name, process, channel, **kwargs)
        else:
            h[0] = self.format(
                h[0],
                channel,
                process,
                rebin_pt=rebin_pt,
                rebin_y=rebin_y,
                normalize=normalize,
                apply_postOp=apply_postOp,
            )
            h[1] = self.format(
                h[1],
                channel,
                process,
                rebin_pt=rebin_pt,
                rebin_y=rebin_y,
                normalize=normalize,
                apply_postOp=apply_postOp,
            )

            h = hh.divideHists(h[0], h[1])
            h = hh.multiplyHists(h, self.ref[channel][process])

            super().add_systematic(h, name, process, channel, **kwargs)

    def add_process(
        self,
        h,
        name,
        channel,
        rebin_pt=True,
        rebin_y=True,
        normalize=True,
        apply_postOp=True,
        **kwargs,
    ):

        h = self.format(
            h,
            channel,
            name,
            rebin_pt=rebin_pt,
            rebin_y=rebin_y,
            normalize=normalize,
            apply_postOp=apply_postOp,
        )
        super().add_process(h, name, channel, **kwargs)

        self.ref[channel][name] = h

    def format(
        self,
        h,
        ch,
        process,
        rebin_pt=True,
        rebin_y=True,
        normalize=True,
        apply_postOp=True,
    ):
        """
        Apply mass, charge selections to input histogram.
        Re-order axes to match the predictions (pt, absY, ...).
        Optionally, re-bin pt and y.
        Returns histogram with changes applied.
        """

        h = copy.deepcopy(h)  # pass by value

        h = self.apply_selections(h, process, ch)

        pt_axis_name = self.get_ptV_axis_name(h)
        absY_axis_name = self.get_absYV_axis_name(h)
        charge_axis_name = self.get_charge_axis_name(h)

        # rename axes to match the reference histogram
        hh.renameAxis(h, pt_axis_name, self.ref[ch]["ptV_name"])
        hh.renameAxis(h, absY_axis_name, self.ref[ch]["absYV_name"])
        if charge_axis_name:
            hh.renameAxis(h, charge_axis_name, self.ref[ch]["chargeV_name"])

        # disable underflow for pt, for compatibility
        h = hh.disableFlow(h, self.ref[ch]["ptV_name"], under=False, over=True)

        # optionally, rebin pt and y, and scale
        if rebin_pt:
            h = hh.rebinHist(h, self.ref[ch]["ptV_name"], self.ref[ch]["ptV_bins"])
        if rebin_y:
            h = hh.rebinHist(h, self.ref[ch]["absYV_name"], self.ref[ch]["absYV_bins"])
        if normalize:
            h *= self.ref[ch]["lumi"] * self.ref[ch]["scale"]

        # re-order remaining axes to match expected output
        remaining_axes = list(h.axes.name)
        remaining_axes.remove(self.ref[ch]["ptV_name"])
        remaining_axes.remove(self.ref[ch]["absYV_name"])
        ordered_axes = [
            self.ref[ch]["ptV_name"],
            self.ref[ch]["absYV_name"],
            *remaining_axes,
        ]
        h = h.project(*ordered_axes)

        # optionally, apply any post operation
        if self.ref[ch]["postOp"] is not None and apply_postOp:
            h = self.ref[ch]["postOp"](h)

        return h

    def get_ptV_axis_name(self, h):
        for name in ["ptVgen", "ptVGen", "qT"]:
            if name in h.axes.name:
                return name
        self.logger.debug(f"Did not find pT axis! Available axes: {h.axes.name}")

    def get_absYV_axis_name(self, h):
        for name in ["absYVgen", "absYVGen", "absY"]:
            if name in h.axes.name:
                return name
        self.logger.debug(f"Did not find absY axis! Available axes: {h.axes.name}")

    def get_ptLep_axis_name(self, h):
        for name in ["ptGen"]:
            if name in h.axes.name:
                return name
        self.logger.debug(f"Did not find pT axis! Available axes: {h.axes.name}")

    def get_absEtaLep_axis_name(self, h):
        for name in ["absEtaGen"]:
            if name in h.axes.name:
                return name
        self.logger.debug(f"Did not find absY axis! Available axes: {h.axes.name}")

    def get_charge_axis_name(self, h):
        for name in ["chargeVgen", "charge", "qGen"]:
            if name in h.axes.name:
                return name
        self.logger.debug(f"Did not find charge axis! Available axes: {h.axes.name}")

    def get_mass_axis_name(self, h):
        for name in ["massVgen", "Q"]:
            if name in h.axes.name:
                return name
        self.logger.debug(f"Did not find mass axis! Available axes: {h.axes.name}")

    def apply_selections(self, h, process, channel):

        if process == "Zmumu":

            # mass selection
            mass_axis_name = self.get_mass_axis_name(h)
            if mass_axis_name:
                h = h[{mass_axis_name: hist.sum}]

            # charge selection
            charge_axis_name = self.get_charge_axis_name(h)
            if charge_axis_name:
                h = h[{charge_axis_name: 0.0j}]

            # select sigmaUL if needed
            if channel == "chSigmaUL" and "helicity" in h.axes.name:
                h = h[{"helicity": -1.0j}]

        elif process == "Wmunu":

            # mass selection
            mass_axis_name = self.get_mass_axis_name(h)
            if mass_axis_name:
                h = h[{mass_axis_name: hist.sum}]

        return h


def calculate_ais_from_helicities_hist(h_hels):
    """
    Calculate A_i histogram from helicities histogram.
    Assuming as input a histogram with the helicities,
    with the first element being the sigma_UL.
    Returns a the A_i's histogram.
    """

    h_hels = copy.deepcopy(h_hels)  # pass by value

    # re-order histogram axes such that the last one is the helicity
    if h_hels.axes[-1].name != "helicity":
        axes = list(h_hels.axes.name)
        axes.remove("helicity")
        axes.append("helicity")
        h_hels = h_hels.project(*axes)

    vals = copy.copy(h_hels.values())
    vars = copy.copy(
        h_hels.variances()
    )  # need to copy in case the storage is hist.storage.Double

    # A_i = sigma_i/sigma_UL
    num = vals[..., 1:]
    den = vals[..., 0][..., np.newaxis]
    vals[..., 1:] = np.where(den == 0, 0.0, num / den)

    # treat these as poisson uncorrelated rv's
    if np.any(vars):

        vars[..., 1:] = np.where(
            vals[..., 0][..., np.newaxis] == 0,
            0.0,
            den**-2 * vars[..., 1:] + num**2 * den**-4 * vars[..., 0][..., np.newaxis],
        )

    h_out = hist.Hist(
        *h_hels.axes[:-1],
        hist.axis.Integer(0, 8, name="ais"),
        storage=hist.storage.Weight() if np.any(vars) else hist.storage.Double(),
    )
    h_out.values()[...] = vals[..., 1:]
    if np.any(vars):
        h_out.variances()[...] = vars[..., 1:]

    return h_out


def convert_WFull_to_LepFiducial(h_W_lep_fiducial, h_W_lep_inclusive):

    def _convert_WFull_to_LepFiducial(h):

        if "qT" in h.axes.name:
            hh.renameAxis(h, "qT", "ptVgen")
            hh.renameAxis(h, "absY", "absYVgen")
            hh.renameAxis(h, "charge", "chargeVgen")

        if h.axes["absYVgen"].traits.underflow:  # not needed, messes up the dimensions
            h = hh.disableFlow(h, "absYVgen", under=False, over=True)
        if h.axes["ptVgen"].traits.underflow:  # not needed, messes up the dimensions
            h = hh.disableFlow(h, "ptVgen", under=False, over=True)

        correction = hh.divideHists(h, h_W_lep_inclusive)
        final_correction = copy.deepcopy(correction)
        final_correction.values(flow=True)[...] = np.ones_like(
            h_W_lep_inclusive.values(flow=True)
        )
        final_correction.values()[...] = correction.values()
        correction = final_correction

        out = hh.multiplyHists(h_W_lep_fiducial, correction)
        out = out.project("absEtaGen", "ptGen", "qGen")

        return out

    return _convert_WFull_to_LepFiducial


analysis_label = Datagroups.analysisLabel(os.path.basename(__file__))
parser, initargs = parsing.common_parser(analysis_label)

parser.add_argument(
    "infile",
    type=str,
    help="Input unfolded fit result for the Z distributions.",
)
parser.add_argument(
    "--fitresultModel", type=str, default="Select helicitySig:slice(0,1)"
)
parser.add_argument(
    "--predGenerator",
    type=str,
    default=f"scetlib_dyturbo",
    help="Generator used for the sigmaUL predictions. "
    "Expect the prediction file to be named as <generator>CorrZ.pkl.lz4. "
    "Expect the prediction histogram to be named as <generator>_hist",
)
parser.add_argument(
    "--constrainAlphaS", action="store_true", help="Constrain alpha_S in the fit."
)
parser.add_argument(
    "--constrainMW", action="store_true", help="Constrain mW in the fit."
)
parser.add_argument("--noFitSigmaUL", action="store_true", help="Don't fit sigmaUL.")
parser.add_argument(
    "--fitAngularCoeffs",
    action="store_true",
    default=False,
    help="Fit the angular coefficients."
    "Predictions of the Ai's from --predAiFile."
    "Predictions of the sigma_UL from infile.",
)
parser.add_argument(
    "--predAiFile",
    type=str,
    default=f"{common.data_dir}/TheoryCorrections/w_z_helicity_xsecs_scetlib_dyturboCorr_maxFiles_m1_alphaSunfoldingBinning_helicity.hdf5",
    help="Gen file used for the Ai predictions."
    "Will be stitched with the --predGenerator file.",
)
parser.add_argument(
    "-W",
    "--fitW",
    action="store_true",
    help="Include W in the fit."
    "Will use --predWFile for predictions and --infileW for the unfolded distribution.",
)
parser.add_argument("--outname", default="carrot", help="output file name")
parser.add_argument(
    "--sparse",
    default=False,
    action="store_true",
    help="Make sparse tensor",
)
parser.add_argument(
    "--systematicType",
    choices=["log_normal", "normal"],
    default="normal",
    help="probability density for systematic variations",
)

args = parser.parse_args()
logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

# Build tensor
writer = AlphaSTheoryFitTW(
    sparse=args.sparse,
    systematic_type="normal" if args.fitAngularCoeffs else args.systematicType,
    allow_negative_expectation=args.fitAngularCoeffs,
)

# load in data histograms and covariance matrix

fitresult, meta = rabbit.io_tools.get_fitresult(args.infile, result="asimov", meta=True)

# covariance across any number of physics models
# h_data_cov = fitresult["physics_models"]["CompositeModel"][
#     "hist_postfit_inclusive_cov"
# ].get()
h_data_cov = fitresult["physics_models"]["Select"]["hist_postfit_inclusive_cov"].get()
writer.add_data_covariance(h_data_cov)  # N.B: run fit with --externalCovariance

# read and initialize sigmaUL channel
if not args.noFitSigmaUL:
    h_data = fitresult["physics_models"]["CompositeModel"]["channels"][
        "Select helicitySig:slice(0,1)_ch0_masked"
    ]["hist_postfit_inclusive"].get()[:, :, 0]
    writer.add_channel(h_data.axes, "chSigmaUL")
    writer.add_data(h_data, "chSigmaUL")
    writer.set_reference("chSigmaUL", h_data)

# if set, read and initialize Ai's channel
if args.fitAngularCoeffs:
    h_data_ai = fitresult["physics_models"]["CompositeModel"]["channels"][
        "AngularCoefficients ch0_masked ptVGen:rebin(0,3,6,9,12,16,20,24,28,33,44)_ch0_masked"
    ]["hist_postfit_inclusive"].get()
    writer.add_channel(h_data_ai.axes, "chAis")
    writer.add_data(h_data_ai, "chAis")
    writer.set_reference("chAis", h_data_ai, postOp=calculate_ais_from_helicities_hist)

# if set, read and initialize W lepton channel
if args.fitW:
    # h_data_prefsrLep = fitresult["physics_models"]["CompositeModel"]["channels"][
    #     "Select_ch1_masked"
    # ]["hist_postfit_inclusive"].get()
    h_data_prefsrLep = fitresult["physics_models"]["Select"]["channels"]["ch0_masked"][
        "hist_postfit_inclusive"
    ].get()
    writer.add_channel(h_data_prefsrLep.axes, "chW")
    writer.add_data(h_data_prefsrLep, "chW")

# add backgrounds

#  Z->mumu, from your favorite generator
if not args.noFitSigmaUL:
    h_sig_sigmaUL = theory_corrections.load_corr_hist(
        f"{common.data_dir}/TheoryCorrections/{args.predGenerator}CorrZ.pkl.lz4",
        "Z",
        f"{args.predGenerator}_hist",
    )
    h_sig_sigmaUL = h_sig_sigmaUL[{"vars": "pdf0"}]  # select baseline variation
    writer.add_process(h_sig_sigmaUL, "Zmumu", "chSigmaUL", signal=False)

# if set, load in the helicity cross sections predictions from MINNLO
if args.fitAngularCoeffs:
    with h5py.File(args.predAiFile, "r") as ff:
        inputs = input_tools.load_results_h5py(ff)
        h_sig_hels = inputs["Z"]
    h_sig_hels = h_sig_hels[
        {"muRfact": 1.0j, "muFfact": 1.0j}
    ]  # select baseline variation
    writer.add_process(h_sig_hels, "Zmumu", "chAis", signal=False)

# if set, load in the W lepton distributions
if args.fitW:

    # MiNNLO(W, lep; fiducial lep)
    predWFiducialFile = "/ceph/submit/data/group/cms/store/user/lavezzo/alphaS//250710_gen_histmaker/w_z_gen_dists.cash_maxFiles_m1_MiNNLO_lepInclusive.hdf5"
    # TODO move this in wremnants-data? possibly only extract the hists we need
    with h5py.File(predWFiducialFile, "r") as h5file:

        lumi = 16800
        inputs = input_tools.load_results_h5py(h5file)

        weight_sum = inputs["WplusmunuPostVFP"]["weight_sum"]
        xsec = inputs["WplusmunuPostVFP"]["dataset"]["xsec"]
        h_Wp_lep_fiducial = inputs["WplusmunuPostVFP"]["output"][
            "nominal_gen_prefsrlep"
        ].get()
        h_Wp_lep_fiducial *= xsec * lumi / weight_sum

        weight_sum = inputs["WminusmunuPostVFP"]["weight_sum"]
        xsec = inputs["WminusmunuPostVFP"]["dataset"]["xsec"]
        h_Wm_lep_fiducial = inputs["WminusmunuPostVFP"]["output"][
            "nominal_gen_prefsrlep"
        ].get()
        h_Wm_lep_fiducial *= xsec * lumi / weight_sum

        # important to set flow=True since we are inclusive in W gen
        h_W_lep_fiducial = hh.addHists(h_Wp_lep_fiducial, h_Wm_lep_fiducial, flow=True)
        h_W_lep_fiducial = hh.disableFlow(
            h_W_lep_fiducial, "ptVgen", under=False, over=True
        )
        h_W_lep_fiducial = hh.disableFlow(
            h_W_lep_fiducial, "absYVgen", under=False, over=True
        )
        h_W_lep_fiducial = hh.disableFlow(
            h_W_lep_fiducial, "absEtaGen", under=False, over=True
        )
        h_W_lep_fiducial = h_W_lep_fiducial.project(
            "ptVgen", "absYVgen", "chargeVgen", "absEtaGen", "ptGen", "qGen"
        )
        h_W_lep_inclusive = h_W_lep_fiducial.project("ptVgen", "absYVgen", "chargeVgen")

        h_W_lep_fiducial = hh.rebinHist(
            h_W_lep_fiducial, "ptGen", h_data_prefsrLep.axes["ptGen"].edges, flow=True
        )
        h_W_lep_fiducial = hh.rebinHist(
            h_W_lep_fiducial,
            "absEtaGen",
            h_data_prefsrLep.axes["absEtaGen"].edges,
            flow=True,
        )

    # <your favorite generator>(W; full)
    h_pred_W_full = theory_corrections.load_corr_hist(
        f"{common.data_dir}/TheoryCorrections/{args.predGenerator}CorrW.pkl.lz4",
        "W",
        f"{args.predGenerator}_hist",
    )
    h_pred_W_full = h_pred_W_full[{"vars": "pdf0"}].project("qT", "absY", "charge")
    h_pred_W_full = hh.disableFlow(h_pred_W_full, "qT", under=False, over=True)
    h_pred_W_full = hh.disableFlow(h_pred_W_full, "absY", under=False, over=True)

    writer.set_reference(
        "chW",
        h_pred_W_full,
        postOp=convert_WFull_to_LepFiducial(h_W_lep_fiducial, h_W_lep_inclusive),
    )
    writer.add_process(h_pred_W_full, "Wmunu", "chW", signal=False)

# add systematic variations

bosons = []
if not args.noFitSigmaUL:
    bosons.append("Z")
if args.fitW:
    bosons.append("W")

# alphaS variation
logger.info("Now at variations for AlphaS")
if args.predGenerator == "scetlib_dyturbo":
    alphas_vars = theory_corrections.load_corr_helpers(
        bosons,
        ["scetlib_dyturboCT18Z_pdfas"],
        make_tensor=False,
        minnlo_ratio=False,
    )

    if not args.noFitSigmaUL:
        writer.add_systematic(
            [
                alphas_vars["Z"]["scetlib_dyturboCT18Z_pdfas"][{"vars": 2}],
                alphas_vars["Z"]["scetlib_dyturboCT18Z_pdfas"][{"vars": 1}],
            ],
            "pdfAlphaS",
            "Zmumu",
            "chSigmaUL",
            noi=not args.constrainAlphaS,
            constrained=args.constrainAlphaS,
            symmetrize="average",
            kfactor=1.5 / 2.0,
            groups=(
                ["pdfCT18Z", "pdfCT18ZAlphaS", "theory"]
                if args.constrainAlphaS
                else ["pdfCT18Z"]
            ),
        )

    # alphaS variations for W come from same as Z
    if args.fitW:
        writer.add_systematic(
            [
                alphas_vars["W"]["scetlib_dyturboCT18Z_pdfas"][{"vars": 2}],
                alphas_vars["W"]["scetlib_dyturboCT18Z_pdfas"][{"vars": 1}],
            ],
            "pdfAlphaS",
            "Wmunu",
            "chW",
            noi=not args.constrainAlphaS,
            constrained=args.constrainAlphaS,
            symmetrize="average",
            kfactor=1.5 / 2.0,
            groups=(
                ["pdfCT18Z", "pdfCT18ZAlphaS", "theory"]
                if args.constrainAlphaS
                else ["pdfCT18Z"]
            ),
        )

    # Ai's alphaS predictions come only from MiNNLO
    if args.fitAngularCoeffs:
        with h5py.File(
            args.predAiFile.replace("w_z_helicity_xsecs", "w_z_gen_dists"), "r"
        ) as ff:
            inputs = input_tools.load_results_h5py(ff)
            alpha_vars_hels = inputs["ZmumuPostVFP"]["output"][
                "nominal_gen_helicity_nominal_gen_pdfCT18ZalphaS002"
            ].get()
        writer.add_systematic(
            [
                alpha_vars_hels[{"alphasVar": "as0120"}],
                alpha_vars_hels[{"alphasVar": "as0116"}],
            ],
            "pdfAlphaS",
            "Zmumu",
            "chAis",
            noi=not args.constrainAlphaS,
            constrained=args.constrainAlphaS,
            symmetrize="average",
            kfactor=1.5 / 2.0,
            groups=(
                ["pdfCT18Z", "pdfCT18ZAlphaS", "theory"]
                if args.constrainAlphaS
                else ["pdfCT18Z"]
            ),
        )

else:
    raise Exception("No valid configuration found for alphaS variation.")

# # mW variations
# if args.fitW:

#     logger.info("Now at mW variations")

#     # mW variations come from MiNNLO
#     with h5py.File(
#         "/ceph/submit/data/group/cms/store/user/lavezzo/alphaS/250718_mW_variations/prefsr_massWeightW_MiNNLO.hdf5",
#         "r",
#     ) as ff:
#         inputs = input_tools.load_results_h5py(ff)
#         mass_vars_Wp = inputs["WplusmunuPostVFP"]["prefsr_massWeightW"].get()
#         mass_vars_Wm = inputs["WminusmunuPostVFP"]["prefsr_massWeightW"].get()

#     mass_vars_W = hh.addHists(mass_vars_Wp, mass_vars_Wm)
#     mass_vars_W = mass_vars_W.project(
#         "absEtaGen", "ptGen", "qGen", "massShift"
#     )  # re-order
#     mass_vars_W = hh.rebinHist(
#         mass_vars_W, "ptGen", h_data_prefsrLep.axes["ptGen"].edges
#     )
#     mass_vars_W = hh.rebinHist(
#         mass_vars_W, "absEtaGen", h_data_prefsrLep.axes["absEtaGen"].edges
#     )
#     mass_vars_W_nom = mass_vars_W[{"massShift": "massShiftW0MeV"}]
#     mass_vars_W_up = mass_vars_W[{"massShift": "massShiftW100MeVUp"}]
#     mass_vars_W_down = mass_vars_W[{"massShift": "massShiftW100MeVDown"}]

#     # calculate variations on the leptons
#     # (MiNNLO_var / MiNNLO_nom) * nominal
#     # where nominal is MiNNLO(W, lep; fiducial) * scetlib(W; full) / MiNNLO(W; full)
#     ratio_MiNNLO_up = hh.divideHists(mass_vars_W_up, mass_vars_W_nom)
#     ratio_MiNNLO_down = hh.divideHists(mass_vars_W_down, mass_vars_W_nom)
#     writer.ref["chW"]["Wmunu"] = hh.disableFlow(
#         writer.ref["chW"]["Wmunu"], "absEtaGen", under=False, over=True
#     )
#     mass_var_up = hh.multiplyHists(ratio_MiNNLO_up, writer.ref["chW"]["Wmunu"])
#     mass_var_down = hh.multiplyHists(ratio_MiNNLO_down, writer.ref["chW"]["Wmunu"])

#     writer.add_systematic(
#         [mass_var_up, mass_var_down],
#         "massShiftW",
#         "Wmunu",
#         "chW",
#         mirror=False,
#         noi=not args.constrainMW,
#         constrained=args.constrainMW,
#         format=False,
#         groups=["ZmassAndWidth", "theory"]
#     )

# logger.info(f"Now at variations from {args.predGenerator}")
# generator_vars = theory_corrections.load_corr_helpers(
#     bosons,
#     [
#         args.predGenerator,
#         f"{args.predGenerator}MSHT20mcrange",
#         f"{args.predGenerator}MSHT20mbrange",
#     ],
#     make_tensor=False,
#     minnlo_ratio=False,
# )
# for proc in generator_vars.keys():  # loop over processes

#     if proc == "Z":
#         proc_name = "Zmumu"
#         ch_name = "chSigmaUL"
#     elif proc == "W":
#         proc_name = "Wmunu"
#         ch_name = "chW"

#     h = generator_vars[proc][args.predGenerator]

#     # correlated NP uncertainties
#     corr_NP_uncs = [
#         ["Lambda20.25", "Lambda2-0.25", "chargeVgenNP0scetlibNPZLambda2"],
#         ["Lambda4.16", "Lambda4.01", "chargeVgenNP0scetlibNPZLambda4"],
#         [
#             "Delta_Lambda20.02",
#             "Delta_Lambda2-0.02",
#             "chargeVgenNP0scetlibNPZDelta_Lambda2",
#         ],
#     ]
#     for var in corr_NP_uncs:
#         writer.add_systematic(
#             [h[{"vars": var[0]}], h[{"vars": var[1]}]],
#             var[2],
#             proc_name,
#             ch_name,
#             symmetrize="average",
#             groups=["resumNonpert", "resum", "pTModeling", "theory"],
#         )

#     # gamma NP uncertainties
#     gamma_NP_uncs = [
#         ["omega_nu0.5", "c_nu-0.1-omega_nu0.5", "scetlibNPgamma"],
#     ]
#     for var in gamma_NP_uncs:
#         writer.add_systematic(
#             [h[{"vars": var[0]}], h[{"vars": var[1]}]],
#             var[2],
#             proc_name,
#             ch_name,
#             symmetrize="average",
#             groups=["resumTNP", "resum", "pTModeling", "theory"],
#         )

#     # TNP
#     TNP_uncs = [
#         ["gamma_cusp1.", "gamma_cusp-1."],
#         ["gamma_mu_q1.", "gamma_mu_q-1."],
#         ["gamma_nu1.", "gamma_nu-1."],
#         ["h_qqV1.", "h_qqV-1."],
#         ["s1.", "s-1."],
#         ["b_qqV0.5", "b_qqV-0.5"],
#         ["b_qqbarV0.5", "b_qqbarV-0.5"],
#         ["b_qqS0.5", "b_qqS-0.5"],
#         ["b_qqDS0.5", "b_qqDS-0.5"],
#         ["b_qg0.5", "b_qg-0.5"],
#     ]
#     for var in TNP_uncs:
#         var_name = "resumTNP_" + var[1].split("-")[0]
#         writer.add_systematic(
#             [h[{"vars": var[0]}], h[{"vars": var[1]}]],
#             var_name,
#             proc_name,
#             ch_name,
#             symmetrize="average",
#             groups=["resumTNP", "resum", "pTModeling", "theory"],
#         )

#     # transition FO scale uncertainties
#     transition_FO_uncs = [
#         [
#             "transition_points0.2_0.75_1.0",
#             "transition_points0.2_0.35_1.0",
#             "resumTransitionZ",
#         ],
#         [
#             "renorm_scale_pt20_envelope_Up",
#             "renorm_scale_pt20_envelope_Down",
#             "resumFOScaleZ",
#         ],
#     ]
#     for var in transition_FO_uncs:
#         writer.add_systematic(
#             [h[{"vars": var[0]}], h[{"vars": var[1]}]],
#             var[2],
#             proc_name,
#             ch_name,
#             symmetrize="quadratic",
#             groups=["resumTransitionFOScale", "resum", "pTModeling", "theory"],
#         )

#     # mass quark effects
#     h = generator_vars[proc][f"{args.predGenerator}MSHT20mbrange"]
#     writer.add_scale_systematic(
#         [h[{"vars": -1}], h[{"vars": 1}], h[{"vars": 0}]],
#         "pdfMSHT20mbrange",
#         proc_name,
#         ch_name,
#         symmetrize="quadratic",
#         groups=["bcQuarkMass", "pTModeling", "theory"],
#     )
#     h = generator_vars[proc][f"{args.predGenerator}MSHT20mcrange"]
#     writer.add_scale_systematic(
#         [h[{"vars": -1}], h[{"vars": 1}], h[{"vars": 0}]],
#         "pdfMSHT20mcrange",
#         proc_name,
#         ch_name,
#         symmetrize="quadratic",
#         groups=["bcQuarkMass", "pTModeling", "theory"],
#     )

# # PDF uncertainties
# logger.info("Now at PDF variations")
# if args.fitAngularCoeffs:
#     # for Ai's, we have MINNLO, so use it for sigmaUL + Ai's to be consistent

#     # TODO fix this at some point
#     with h5py.File(
#         # args.predAiFile.replace("w_z_helicity_xsecs", "w_z_gen_dists"), "r"
#         "/ceph/submit/data/group/cms/store/user/lavezzo/alphaS//250627_angularCoefficients/w_z_gen_dists.cash_scetlib_dyturboCorr_maxFiles_1000_alphaSunfoldingBinning_helicity_WZ.hdf5",
#         "r",
#     ) as ff:
#         inputs = input_tools.load_results_h5py(ff)
#         pdf_vars = inputs["ZmumuPostVFP"]["output"][
#             "nominal_gen_helicity_pdfCT18Z"
#         ].get()
#         pdf_vars_Wp = inputs["WplusmunuPostVFP"]["output"][
#             "nominal_gen_helicity_pdfCT18Z"
#         ].get()
#         pdf_vars_Wm = inputs["WminusmunuPostVFP"]["output"][
#             "nominal_gen_helicity_pdfCT18Z"
#         ].get()
#         pdf_vars_W = hh.addHists(pdf_vars_Wp, pdf_vars_Wm)

#     for ivar in range(1, len(pdf_vars.axes[-1]), 2):

#         # sigmaUL
#         if not args.noFitSigmaUL:
#             writer.add_scale_systematic(
#                 [
#                     pdf_vars[{"pdfVar": ivar + 1}],
#                     pdf_vars[{"pdfVar": ivar}],
#                     pdf_vars[{"pdfVar": "pdf0CT18Z"}],
#                 ],
#                 f"pdf{int((ivar+1)/2)}CT18Z",
#                 "Zmumu",
#                 "chSigmaUL",
#                 symmetrize="quadratic",
#                 kfactor=1 / 1.645,
#                 groups=["pdfCT18Z", f"pdfCT18ZNoAlphaS", "theory"],
#             )

#         # Ai's
#         writer.add_scale_systematic(
#             [
#                 pdf_vars[{"pdfVar": ivar + 1}],
#                 pdf_vars[{"pdfVar": ivar}],
#                 pdf_vars[{"pdfVar": "pdf0CT18Z"}],
#             ],
#             f"pdf{int((ivar+1)/2)}CT18Z",
#             "Zmumu",
#             "chAis",
#             symmetrize="quadratic",
#             kfactor=1 / 1.645,
#             groups=["pdfCT18Z", f"pdfCT18ZNoAlphaS", "theory"],
#         )

#         if args.fitW:
#             # TODO I'm 90% sure this is wrong
#             raise Exception()
#             writer.add_systematic(
#                 [
#                     pdf_vars_W[{"helicity": -1j}][{"pdfVar": ivar + 1}],
#                     pdf_vars_W[{"helicity": -1j}][{"pdfVar": ivar}],
#                     pdf_vars_W[{"helicity": -1j}][{"pdfVar": "pdf0CT18Z"}],
#                 ],
#                 f"pdf{int((ivar+1)/2)}CT18Z",
#                 "Wmunu",
#                 "chW",
#                 symmetrize="quadratic",
#                 kfactor=1 / 1.645,
#                 groups=["pdfCT18Z", f"pdfCT18ZNoAlphaS", "theory"],
#             )

# else:

#     # for sigmaUL only, scetlib+dyturbo has the latest & greatest PDFs
#     corr_helpers = theory_corrections.load_corr_helpers(
#         bosons,
#         ["scetlib_dyturboCT18ZVars"],
#         make_tensor=False,
#         minnlo_ratio=False,
#     )

#     if not args.noFitSigmaUL:
#         h = corr_helpers["Z"]["scetlib_dyturboCT18ZVars"]
#         for ivar in range(1, len(h.axes[-1]), 2):
#             writer.add_systematic(
#                 [h[{"vars": ivar + 1}], h[{"vars": ivar}]],
#                 f"pdf{int((ivar+1)/2)}CT18Z",
#                 "Zmumu",
#                 "chSigmaUL",
#                 symmetrize="quadratic",
#                 kfactor=1 / 1.645,
#                 groups=["pdfCT18Z", f"pdfCT18ZNoAlphaS", "theory"],
#             )

#     if args.fitW:
#         h = corr_helpers["W"]["scetlib_dyturboCT18ZVars"]
#         for ivar in range(1, len(h.axes[-1]), 2):
#             writer.add_systematic(
#                 [h[{"vars": ivar + 1}], h[{"vars": ivar}]],
#                 f"pdf{int((ivar+1)/2)}CT18Z",
#                 "Wmunu",
#                 "chW",
#                 symmetrize="quadratic",
#                 kfactor=1 / 1.645,
#                 groups=["pdfCT18Z", f"pdfCT18ZNoAlphaS", "theory"],
#             )

# Ai's only uncertainties
if args.fitAngularCoeffs:

    qcd_helper = theory_corrections.make_qcd_uncertainty_helper_by_helicity(
        is_z=True,
        filename=args.predAiFile,
        rebin_ptVgen=h_data_ai.axes["ptVGen"].edges.tolist(),
        rebin_absYVgen=h_data_ai.axes["absYVGen"].edges.tolist(),
        rebin_massVgen=True,
        return_tensor=False,
    )

    # pythia showering uncertainties
    logger.info("Now at pythia_shower_kt")
    pythia_shower_kt = qcd_helper[{"vars": "pythia_shower_kt"}]
    writer.add_scale_systematic(
        [pythia_shower_kt[{"corr": 1}], pythia_shower_kt[{"corr": 0}]],
        "pythia_shower_kt",
        "Zmumu",
        "chAis",
        mirror=True,
        groups=["helicity_shower_kt", "angularCoeffs", "theory"],
    )

    # QCD scales
    logger.info("Now at QCD scales")

    # prepare fine binning hists
    fine_pt_binning = qcd_helper.axes["ptVgen"].edges
    nptfine = len(fine_pt_binning) - 1
    scale_inclusive = np.sqrt((nptfine - 1) / nptfine)

    for hel in range(0, 7 + 1):  # no correction on sigma_UL

        # fine binning
        qcd_scales_hel_up = qcd_helper[{"vars": f"helicity_{hel}_Up"}].project(
            "ptVgen", "absYVgen", "helicity"
        )
        qcd_scales_hel_down = qcd_helper[{"vars": f"helicity_{hel}_Down"}].project(
            "ptVgen", "absYVgen", "helicity"
        )
        qcd_scales_hel_nominal = qcd_helper[{"vars": f"nominal"}].project(
            "ptVgen", "absYVgen", "helicity"
        )

        for bin in range(len(fine_pt_binning) - 1):

            ptl = fine_pt_binning[bin]
            pth = fine_pt_binning[bin + 1]

            qcd_scales_hel_pt_up = copy.deepcopy(qcd_scales_hel_up)
            qcd_scales_hel_pt_up.values()[...] = qcd_scales_hel_nominal.values()
            qcd_scales_hel_pt_up.values()[bin, ...] = qcd_scales_hel_up[
                {"ptVgen": bin}
            ].values()

            qcd_scales_hel_pt_down = copy.deepcopy(qcd_scales_hel_down)
            qcd_scales_hel_pt_down.values()[...] = qcd_scales_hel_nominal.values()
            qcd_scales_hel_pt_down.values()[bin, ...] = qcd_scales_hel_up[
                {"ptVgen": bin}
            ].values()

            writer.add_systematic(
                [qcd_scales_hel_pt_up, qcd_scales_hel_pt_down],
                f"QCDscaleZfine_Pt{ptl}_{pth}helicity_{hel}",
                "Zmumu",
                "chAis",
                symmetrize="quadratic",
                groups=["QCDScaleZMiNNLO", "QCDscale", "angularCoeffs", "theory"],
            )

        # inclusive
        inclusive_pt_binning = [
            qcd_scales_hel_up.axes["ptVgen"].edges[0],
            qcd_scales_hel_up.axes["ptVgen"].edges[-1],
        ]
        qcd_scales_hel_int_up = hh.rebinHist(
            qcd_scales_hel_up, "ptVgen", inclusive_pt_binning
        )
        qcd_scales_hel_int_down = hh.rebinHist(
            qcd_scales_hel_down, "ptVgen", inclusive_pt_binning
        )

        for bin in range(len(fine_pt_binning) - 1):
            qcd_scales_hel_up.values()[bin, ...] = qcd_scales_hel_int_up.values()
            qcd_scales_hel_down.values()[bin, ...] = qcd_scales_hel_int_down.values()

        writer.add_systematic(
            [qcd_scales_hel_up, qcd_scales_hel_down],
            f"QCDscaleZinclusive_Pt{inclusive_pt_binning[0]}_{inclusive_pt_binning[1]}helicity_{hel}",
            "Zmumu",
            "chAis",
            kfactor=scale_inclusive,
            symmetrize="quadratic",
            groups=["QCDScaleZMiNNLO", "QCDscale", "angularCoeffs", "theory"],
        )


# QCD uncertainties by helicity
if args.fitW:

    qcd_helper = theory_corrections.make_qcd_uncertainty_helper_by_helicity(
        is_z=False,
        filename=f"{common.data_dir}/angularCoefficients/w_z_helicity_xsecs_scetlib_dyturboCorr_maxFiles_m1_unfoldingBinning_W.hdf5",
        rebin_ptVgen=False,  # h_pred_W_full.axes["qT"].edges.tolist(),
        rebin_absYVgen=False,  # h_pred_W_full.axes["absY"].edges.tolist(),
        rebin_massVgen=True,
        return_tensor=False,
    )

    print(qcd_helper)
    base = qcd_helper[{"vars": "pythia_shower_kt"}][{"corr": 0}]
    comp = h_pred_W_full
    comp = hh.rebinHist(comp, "qT", base.axes["ptVgen"].edges)
    comp = hh.rebinHist(comp, "absY", base.axes["absYVgen"].edges)
    print(base.values() / comp.values())

    # pythia showering uncertainties
    logger.info("Now at pythia_shower_kt")
    pythia_shower_kt = qcd_helper[{"vars": "pythia_shower_kt"}]
    writer.add_scale_systematic(
        [pythia_shower_kt[{"corr": 1}], pythia_shower_kt[{"corr": 0}]],
        "pythia_shower_kt",
        "Zmumu",
        "chAis",
        mirror=True,
        groups=["helicity_shower_kt", "angularCoeffs", "theory"],
    )

    # QCD scales
    logger.info("Now at QCD scales")

    # prepare fine binning hists
    fine_pt_binning = qcd_helper.axes["ptVgen"].edges
    nptfine = len(fine_pt_binning) - 1
    scale_inclusive = np.sqrt((nptfine - 1) / nptfine)

    for hel in range(0, 7 + 1):  # no correction on sigma_UL

        # fine binning
        qcd_scales_hel_up = qcd_helper[{"vars": f"helicity_{hel}_Up"}].project(
            "ptVgen", "absYVgen", "helicity"
        )
        qcd_scales_hel_down = qcd_helper[{"vars": f"helicity_{hel}_Down"}].project(
            "ptVgen", "absYVgen", "helicity"
        )
        qcd_scales_hel_nominal = qcd_helper[{"vars": f"nominal"}].project(
            "ptVgen", "absYVgen", "helicity"
        )

        for bin in range(len(fine_pt_binning) - 1):

            ptl = fine_pt_binning[bin]
            pth = fine_pt_binning[bin + 1]

            qcd_scales_hel_pt_up = copy.deepcopy(qcd_scales_hel_up)
            qcd_scales_hel_pt_up.values()[...] = qcd_scales_hel_nominal.values()
            qcd_scales_hel_pt_up.values()[bin, ...] = qcd_scales_hel_up[
                {"ptVgen": bin}
            ].values()

            qcd_scales_hel_pt_down = copy.deepcopy(qcd_scales_hel_down)
            qcd_scales_hel_pt_down.values()[...] = qcd_scales_hel_nominal.values()
            qcd_scales_hel_pt_down.values()[bin, ...] = qcd_scales_hel_up[
                {"ptVgen": bin}
            ].values()

            writer.add_systematic(
                [qcd_scales_hel_pt_up, qcd_scales_hel_pt_down],
                f"QCDscaleZfine_Pt{ptl}_{pth}helicity_{hel}",
                "Zmumu",
                "chAis",
                symmetrize="quadratic",
                groups=["QCDScaleZMiNNLO", "QCDscale", "angularCoeffs", "theory"],
            )

        # inclusive
        inclusive_pt_binning = [
            qcd_scales_hel_up.axes["ptVgen"].edges[0],
            qcd_scales_hel_up.axes["ptVgen"].edges[-1],
        ]
        qcd_scales_hel_int_up = hh.rebinHist(
            qcd_scales_hel_up, "ptVgen", inclusive_pt_binning
        )
        qcd_scales_hel_int_down = hh.rebinHist(
            qcd_scales_hel_down, "ptVgen", inclusive_pt_binning
        )

        for bin in range(len(fine_pt_binning) - 1):
            qcd_scales_hel_up.values()[bin, ...] = qcd_scales_hel_int_up.values()
            qcd_scales_hel_down.values()[bin, ...] = qcd_scales_hel_int_down.values()

        writer.add_systematic(
            [qcd_scales_hel_up, qcd_scales_hel_down],
            f"QCDscaleZinclusive_Pt{inclusive_pt_binning[0]}_{inclusive_pt_binning[1]}helicity_{hel}",
            "Zmumu",
            "chAis",
            kfactor=scale_inclusive,
            symmetrize="quadratic",
            groups=["QCDScaleZMiNNLO", "QCDscale", "angularCoeffs", "theory"],
        )


# write output
directory = args.outfolder
if directory == "":
    directory = "./"
filename = args.outname
if not args.noFitSigmaUL:
    filename += f"_sigmaUL"
if args.fitAngularCoeffs:
    filename += f"_Ais"
if args.fitW:
    filename += f"_W"
if args.postfix:
    filename += f"_{args.postfix}"
writer.write(outfolder=directory, outfilename=filename)
logger.info(f"Written to {os.path.join(directory, filename)}.hdf5")
