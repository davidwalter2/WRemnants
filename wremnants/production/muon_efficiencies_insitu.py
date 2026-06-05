import pickle

import hist
import lz4.frame
import ROOT

import narf
from wremnants.production.muon_efficiencies_smooth import cloneAxis
from wremnants.utilities import common
from wums import logging

logger = logging.child_logger(__name__)

narf.clingutils.Declare('#include "muon_efficiencies_insitu.hpp"')

data_dir = common.data_dir

# In-situ efficiency steps floated as Chebyshev polynomials.
# IDIP   : 1D in pt,           charge-dependent  (NCoeffPt per (eta, q))
# Trigger: 2D in (pt, ut),     charge-dependent  (NCoeffPt*NCoeffUt per (eta,q))
# Iso    : 2D in (pt, ut),     charge-inclusive  (NCoeffPt*NCoeffUt per eta)
# matching the SMP-23-002 measurement variables (pt^3 x ut^2 polynomials).
insitu_eff_steps = ["idip", "trigger", "iso"]

# 48-bin probe-eta decorrelation; probe pt window [26, 60] (34 1-GeV bins).
# ut window [-30, +100] matches the published smoothSF3D_uTm30to100 file.
insitu_n_eta = 48
insitu_n_coeff_pt = 4  # 3rd-order Chebyshev in pt
insitu_n_coeff_ut = 3  # 2nd-order Chebyshev in ut (HLT, Iso only)
insitu_pt_range = (26.0, 60.0)
insitu_ut_range = (-30.0, 100.0)
# variation step folded into the stored response (units of the unconstrained
# nuisance: theta_c = delta * n_c). Small to avoid exp() overflow from the
# large fail-leg derivative of high-efficiency steps.
insitu_delta = 0.01
insitu_effMC_max = 0.999

# group-name prefix per step (used for nuisance grouping in setupRabbit)
insitu_step_group = {
    "idip": "effInsituID",
    "trigger": "effInsituHLT",
    "iso": "effInsituIso",
}


def _build_effMC_hist(h_in, axes_spec, name):
    """Clone an input boost hist into a Weight-storage hist with flow filled
    from the adjacent in-range bin. ``axes_spec`` is a list of
    (axis_name, in_axis, with_flow) tuples giving the desired output axes
    (with flow as requested). The hist is filled positionally from h_in.
    """
    out_axes = []
    flow_dims = []
    for ax_name, ax_in, with_flow in axes_spec:
        if isinstance(ax_in, hist.axis.StrCategory):
            out_axes.append(hist.axis.StrCategory(list(ax_in), name=ax_name))
            flow_dims.append(False)
        else:
            out_axes.append(
                cloneAxis(
                    ax_in,
                    overflow=with_flow,
                    underflow=with_flow,
                    newName=ax_name,
                )
            )
            flow_dims.append(with_flow)
    h = hist.Hist(*out_axes, name=name, storage=hist.storage.Weight())
    h.view(flow=False)["value"] = h_in.values(flow=False)
    # Set under/overflow per axis to the adjacent in-range bin
    view = h.view(flow=True)
    for i_axis, has_flow in enumerate(flow_dims):
        if not has_flow:
            continue
        size_with_flow = out_axes[i_axis].extent
        # underflow -> first in-range
        idx_under = [slice(None)] * view.ndim
        idx_under[i_axis] = 0
        idx_first = [slice(None)] * view.ndim
        idx_first[i_axis] = 1
        view[tuple(idx_under)] = view[tuple(idx_first)]
        # overflow -> last in-range
        idx_over = [slice(None)] * view.ndim
        idx_over[i_axis] = size_with_flow - 1
        idx_last = [slice(None)] * view.ndim
        idx_last[i_axis] = size_with_flow - 2
        view[tuple(idx_over)] = view[tuple(idx_last)]
    return h


def make_muon_insitu_effMC_helper(effMCFile):
    """Load the per-step MC-efficiency hists for the in-situ method.

    ``effMCFile`` is the .pkl.lz4 produced by scripts/corrections/
    make_insitu_effMC.py via wums write_lz4_pkl_output. It holds a dict
        { <process>: { "process": str, "effMC": { "idip": Hist, "trigger":
          Hist, "iso": Hist }, "steps": [...] }, "meta_data": {...},
          "file_meta_data": {...} }
    i.e. the effMC content is nested under a top-level process key alongside
    the provenance metadata. Here:
      - idip   : axes (eta(48), pt(34), charge(2)),                 3D
      - trigger: axes (eta(48), pt(34), charge(2), ut(N_ut)),       4D
      - iso    : axes (eta(48), pt(34), ut(N_ut)),                  3D (charge-incl)

    Returns a dict {step: hist.Hist with flow filled from adjacent bins} so
    the C++ helper does direct per-bin lookups without ever reading an empty
    flow cell.
    """
    with lz4.frame.open(effMCFile, "rb") as fin:
        payload = pickle.load(fin)
    # write_lz4_pkl_output nests the content under a top-level process key,
    # alongside the meta_data / file_meta_data provenance entries.
    content = next(
        v
        for k, v in payload.items()
        if k not in ("meta_data", "file_meta_data", "meta_info")
    )
    effMC_in = content["effMC"]
    assert isinstance(effMC_in, dict) and set(effMC_in.keys()) >= set(
        insitu_eff_steps
    ), f"effMC pkl must be a dict with keys {insitu_eff_steps}, got {list(effMC_in)}"
    logger.info(
        f"Loaded in-situ effMC from {effMCFile} "
        f"(process={content.get('process')}, shapes="
        f"{ {k: effMC_in[k].shape for k in insitu_eff_steps} })"
    )

    h_idip_in = effMC_in["idip"]
    h_trig_in = effMC_in["trigger"]
    h_iso_in = effMC_in["iso"]

    eta_in_idip, pt_in_idip, q_in_idip = h_idip_in.axes
    eta_in_trig, pt_in_trig, q_in_trig, ut_in_trig = h_trig_in.axes
    eta_in_iso, pt_in_iso, ut_in_iso = h_iso_in.axes

    h_idip = _build_effMC_hist(
        h_idip_in,
        [
            ("SF eta", eta_in_idip, True),
            ("SF pt", pt_in_idip, True),
            ("SF charge", q_in_idip, False),
        ],
        name="effMC_insitu_idip",
    )
    h_trig = _build_effMC_hist(
        h_trig_in,
        [
            ("SF eta", eta_in_trig, True),
            ("SF pt", pt_in_trig, True),
            ("SF charge", q_in_trig, False),
            ("SF ut", ut_in_trig, True),
        ],
        name="effMC_insitu_trigger",
    )
    h_iso = _build_effMC_hist(
        h_iso_in,
        [
            ("SF eta", eta_in_iso, True),
            ("SF pt", pt_in_iso, True),
            ("SF ut", ut_in_iso, True),
        ],
        name="effMC_insitu_iso",
    )
    logger.info(
        f"Built in-situ effMC hists: idip={h_idip.shape}, trigger={h_trig.shape}, "
        f"iso={h_iso.shape}"
    )
    return {"idip": h_idip, "trigger": h_trig, "iso": h_iso}


def insitu_parameter_labels(
    n_eta=insitu_n_eta,
    n_coeff_pt=insitu_n_coeff_pt,
    n_coeff_ut=insitu_n_coeff_ut,
):
    """Ordered list of nuisance labels, one per flat tensor index.

    Layout must match muon_efficiencies_insitu.hpp:
      [ idip (q x eta x kPt)
        trigger (q x eta x kPt x mUt)
        iso     (eta x kPt x mUt) ]
    """
    labels = []
    # idip: charge x eta x pt
    for qbit in range(2):
        qtag = "plus" if qbit else "minus"
        for b in range(n_eta):
            for k in range(n_coeff_pt):
                labels.append(f"{insitu_step_group['idip']}_eta{b}_q{qtag}_cPt{k}")
    # trigger: charge x eta x pt x ut
    for qbit in range(2):
        qtag = "plus" if qbit else "minus"
        for b in range(n_eta):
            for k in range(n_coeff_pt):
                for m in range(n_coeff_ut):
                    labels.append(
                        f"{insitu_step_group['trigger']}_eta{b}_q{qtag}"
                        f"_cPt{k}_cUt{m}"
                    )
    # iso: eta x pt x ut (charge-inclusive)
    for b in range(n_eta):
        for k in range(n_coeff_pt):
            for m in range(n_coeff_ut):
                labels.append(f"{insitu_step_group['iso']}_eta{b}_cPt{k}_cUt{m}")
    return labels


def make_muon_insitu_efficiency_helper(
    effMCFile,
    n_coeff_pt=insitu_n_coeff_pt,
    n_coeff_ut=insitu_n_coeff_ut,
    pt_range=insitu_pt_range,
    ut_range=insitu_ut_range,
    delta=insitu_delta,
    effMC_max=insitu_effMC_max,
    single_leg=False,
):
    """Build the RDF helper that emits the in-situ efficiency syst tensor.

    The number of eta bins is inferred from the effMC histograms so it follows
    the histmaker --eta binning automatically. Returns (helper, labels).

    ``single_leg=True`` selects the W single-muon helper (no tag leg); the
    default two-leg helper is used for the Z dilepton tag-and-probe topology.
    """
    effMC = make_muon_insitu_effMC_helper(effMCFile)
    h_idip = effMC["idip"]
    h_trig = effMC["trigger"]
    h_iso = effMC["iso"]
    n_eta = h_idip.axes[0].size
    assert h_trig.axes[0].size == n_eta and h_iso.axes[0].size == n_eta, (
        n_eta,
        h_trig.axes[0].size,
        h_iso.axes[0].size,
    )

    idip_pyroot = narf.hist_to_pyroot_boost(h_idip)
    trig_pyroot = narf.hist_to_pyroot_boost(h_trig)
    iso_pyroot = narf.hist_to_pyroot_boost(h_iso)

    helper_class = (
        ROOT.wrem.muon_insitu_efficiency_helper_singleleg
        if single_leg
        else ROOT.wrem.muon_insitu_efficiency_helper
    )
    helper = helper_class[
        n_eta,
        n_coeff_pt,
        n_coeff_ut,
        type(idip_pyroot),
        type(trig_pyroot),
        type(iso_pyroot),
    ](
        ROOT.std.move(idip_pyroot),
        ROOT.std.move(trig_pyroot),
        ROOT.std.move(iso_pyroot),
        pt_range[0],
        pt_range[1],
        ut_range[0],
        ut_range[1],
        delta,
        effMC_max,
    )

    n_id = n_eta * 2 * n_coeff_pt
    n_hlt = n_eta * 2 * n_coeff_pt * n_coeff_ut
    n_iso_p = n_eta * n_coeff_pt * n_coeff_ut
    n_sf = n_id + n_hlt + n_iso_p
    axis_insitu = hist.axis.Integer(
        0, n_sf, underflow=False, overflow=False, name="insituEffParm"
    )
    helper.tensor_axes = [axis_insitu]

    labels = insitu_parameter_labels(n_eta, n_coeff_pt, n_coeff_ut)
    assert len(labels) == n_sf, (len(labels), n_sf)
    logger.info(
        f"Built in-situ efficiency helper with {n_sf} unconstrained "
        f"Chebyshev coefficients "
        f"(eta={n_eta}, order_pt={n_coeff_pt - 1}, order_ut={n_coeff_ut - 1}; "
        f"idip={n_id}, trigger={n_hlt}, iso={n_iso_p})"
    )
    return helper, labels
