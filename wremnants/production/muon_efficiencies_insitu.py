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

# In-situ efficiency steps floated as Chebyshev polynomials in pt.
# "idip" and "trigger" are charge dependent (see
# common.muonEfficiency_chargeDependentSteps); "iso" is charge inclusive.
insitu_eff_steps = ["idip", "trigger", "iso"]

# Chebyshev order 3 -> 4 coefficients; 48-bin probe-eta decorrelation;
# probe pt window [26, 60] (matches the 34-bin pt axis in mz_dilepton.py).
insitu_n_eta = 48
insitu_n_coeff = 4
insitu_pt_range = (26.0, 60.0)
# variation step folded into the stored response (units of the unconstrained
# nuisance: theta_c = delta * n_c). Small to avoid exp() overflow from the
# large fail-leg derivative of high-efficiency steps; physics is delta
# independent since the coefficients are unconstrained.
insitu_delta = 0.01
insitu_effMC_max = 0.999

# group-name prefix per step (used for nuisance grouping in setupRabbit)
insitu_step_group = {
    "idip": "effInsituID",
    "trigger": "effInsituHLT",
    "iso": "effInsituIso",
}


def make_muon_insitu_effMC_helper(effMCFile):
    """Load the MC-efficiency boost histogram for the in-situ method.

    ``effMCFile`` is the .pkl.lz4 produced by scripts/corrections/
    make_insitu_effMC.py from our own Z->mumu signal MC. It holds a
    hist.Hist with axes [eta(48,[-2.4,2.4]), pt(34,[26,60]),
    charge(Regular(2,-2,2)), insitu_step(StrCategory["idip","trigger","iso"])].

    The eta/pt binning is the Chebyshev decorrelation grid and probe pt
    window by construction, so the C++ helper does a direct per-bin lookup.
    This returns an equivalent hist with eta/pt under/overflow bins filled
    from the adjacent in-range bin so the helper never reads an empty flow
    bin for the fail-leg factor.
    """
    with lz4.frame.open(effMCFile, "rb") as fin:
        payload = pickle.load(fin)
    effMC_in = payload["effMC"]
    logger.info(
        f"Loaded in-situ effMC from {effMCFile} "
        f"(process={payload.get('process')}, shape={effMC_in.shape})"
    )

    ax_eta_in, ax_pt_in, ax_charge_in, ax_step_in = effMC_in.axes
    axis_eta = cloneAxis(ax_eta_in, overflow=True, underflow=True, newName="SF eta")
    axis_pt = cloneAxis(ax_pt_in, overflow=True, underflow=True, newName="SF pt")
    # charge axis matches the convention in muon_efficiencies_smooth
    # (bin 0 = charge -1, bin 1 = charge +1); no flow
    axis_charge = hist.axis.Regular(
        2, -2.0, 2.0, underflow=False, overflow=False, name="SF charge"
    )
    axis_step = hist.axis.StrCategory(list(ax_step_in), name="insitu_step")

    # Weight storage + .at(...).value() in C++ mirrors the proven
    # muon_efficiencies_smooth helper path (storage type affects how the C++
    # helper reads cells).
    effMC = hist.Hist(
        axis_eta,
        axis_pt,
        axis_charge,
        axis_step,
        name="effMC_insitu",
        storage=hist.storage.Weight(),
    )
    # copy in-range content (input axes have no flow)
    effMC.view(flow=False)["value"] = effMC_in.values(flow=False)

    # set eta/pt under/overflow equal to the adjacent in-range bin
    view = effMC.view(flow=True)
    view[0, ...] = view[1, ...]
    view[axis_eta.extent - 1, ...] = view[axis_eta.extent - 2, ...]
    view[:, 0, ...] = view[:, 1, ...]
    view[:, axis_pt.extent - 1, ...] = view[:, axis_pt.extent - 2, ...]

    logger.info(f"Built in-situ effMC hist {effMC.shape} (eta x pt x charge x step)")
    return effMC


def insitu_parameter_labels(n_eta=insitu_n_eta, n_coeff=insitu_n_coeff):
    """Ordered list of nuisance labels, one per flat tensor index.

    Layout must match muon_efficiencies_insitu.hpp:
      [ idip (charge split) | trigger | iso ], each (eta, coeff) C-ordered.
    Names: effInsitu{ID,HLT,Iso}_eta<b>[_q{plus,minus}]_c<k>.
    """
    labels = []
    # idip: charge split (qbit 0 = minus, 1 = plus)
    for qbit in range(2):
        qtag = "plus" if qbit else "minus"
        for b in range(n_eta):
            for k in range(n_coeff):
                labels.append(f"{insitu_step_group['idip']}_eta{b}_q{qtag}_c{k}")
    # trigger: charge inclusive
    for b in range(n_eta):
        for k in range(n_coeff):
            labels.append(f"{insitu_step_group['trigger']}_eta{b}_c{k}")
    # iso: charge inclusive
    for b in range(n_eta):
        for k in range(n_coeff):
            labels.append(f"{insitu_step_group['iso']}_eta{b}_c{k}")
    return labels


def make_muon_insitu_efficiency_helper(
    effMCFile,
    n_coeff=insitu_n_coeff,
    pt_range=insitu_pt_range,
    delta=insitu_delta,
    effMC_max=insitu_effMC_max,
):
    """Build the RDF helper that emits the in-situ efficiency syst tensor.

    The number of eta bins is inferred from the effMC histogram so it follows
    the histmaker --eta binning automatically. Returns (helper, labels) where
    ``helper`` has ``tensor_axes`` set to a single flat ``insituEffParm`` axis
    and ``labels`` is the ordered list of per-index nuisance names (see
    insitu_parameter_labels).
    """
    effMC = make_muon_insitu_effMC_helper(effMCFile)
    n_eta = effMC.axes[0].size  # eta axis, == histmaker --eta binning
    effMC_pyroot = narf.hist_to_pyroot_boost(effMC)

    helper = ROOT.wrem.muon_insitu_efficiency_helper[
        n_eta, n_coeff, type(effMC_pyroot)
    ](ROOT.std.move(effMC_pyroot), pt_range[0], pt_range[1], delta, effMC_max)

    n_sf = n_eta * 2 * n_coeff + n_eta * n_coeff + n_eta * n_coeff
    axis_insitu = hist.axis.Integer(
        0, n_sf, underflow=False, overflow=False, name="insituEffParm"
    )
    helper.tensor_axes = [axis_insitu]

    labels = insitu_parameter_labels(n_eta, n_coeff)
    assert len(labels) == n_sf, (len(labels), n_sf)
    logger.info(
        f"Built in-situ efficiency helper with {n_sf} unconstrained "
        f"Chebyshev coefficients (eta={n_eta}, order={n_coeff - 1})"
    )
    return helper, labels
