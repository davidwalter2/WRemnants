import os
import pickle

import hist
import lz4.frame
import numpy as np
import ROOT

import narf
import narf.clingutils
from wremnants.utilities import common
from wums import logging

logger = logging.child_logger(__name__)

# Declares wrem::cvh_muon_in_bad_module / wrem::cvh_muon_sf, the two possible
# treatments of the CVH-refit efficiency holes of the badly aligned modules
# (TIB-L2 detId 369141860, fixed in WMass/cmssw eb96caef, and TID-2 detId
# 402666798, still uncorrected). See the header for the measurement.
narf.clingutils.Declare('#include "muon_efficiencies_cvh.hpp"')

data_dir = common.data_dir


def apply_bad_module_veto(
    df, ptCut=15.0, etaCut=2.4, ptMargin=1.0, filter_name="badCvhModuleVeto"
):
    """Reject events with a muon crossing one of the pathological modules.

    Applied to DATA AND MC alike: the affected (eta, phi') rectangles are simply
    removed from the measurement, which makes the data-only CVH refit
    inefficiency irrelevant instead of correcting for it.

    The decision uses the uncorrected muon kinematics (Muon_pt/eta/phi/charge)
    rather than the corrected ones, because Muon_corrected* is an alias of the
    CVH refit and is undefined (charge = -99) exactly for the muons that fall in
    the holes -- in data those muons are absent from 'vetoMuons' altogether,
    so a veto built on them would fire in MC only.

    For the same reason the muon selection here is deliberately looser than the
    veto muon definition (only Muon_looseId and the kinematic acceptance, no
    global/tracker quality since that too is evaluated on corrected quantities):
    it has to be a superset of 'vetoMuons', so that every event whose selection
    could be changed by a failed refit is removed. ptMargin lowers the pt
    threshold to absorb the (sub-percent) difference between corrected and
    uncorrected pt at the veto threshold.
    """
    df = df.Define(
        "Muon_inBadCvhModule",
        "wrem::cvh_muons_in_bad_module(Muon_eta, Muon_phi, Muon_charge, Muon_pt)",
    )
    df = df.Filter(
        f"!ROOT::VecOps::Any(Muon_inBadCvhModule && Muon_looseId"
        f" && Muon_pt > {ptCut - ptMargin} && abs(Muon_eta) < {etaCut})",
        filter_name,
    )
    return df


def define_cvh_weight(df, muons, out_col="weight_cvhSF"):
    """Define the per-event CVH efficiency weight = product of the per-muon SF
    over all reconstructed muons in the final state (1 for W, 2 for Z).

    Alternative to apply_bad_module_veto, not to be combined with it.

    muons: iterable of (eta_expr, phi_expr, charge_expr, pt_expr) tuples, each a
    column name or an inline RDF expression, e.g.
    ("trigMuons_eta0", "trigMuons_phi0", "trigMuons_charge0", "trigMuons_pt0").
    charge and pt are needed to undo the track bending, i.e. to evaluate the SF
    map in the module azimuth phi' = phi - q*C/pt (see muon_efficiencies_cvh.hpp).
    Call on MC only; multiply the result into the nominal weight.
    """
    expr = "*".join(
        f"wrem::cvh_muon_sf({eta}, {phi}, {charge}, {pt})"
        for eta, phi, charge, pt in muons
    )
    df = df.Define(out_col, expr)
    return df, out_col


def _cubic_spline(x, y, xnew):
    """Not-a-knot cubic spline through (x, y), evaluated at xnew.

    Same spline as scipy's CubicSpline / RegularGridInterpolator(method="cubic")
    used in smoothVetoSF.py, reimplemented because scipy.linalg cannot be called
    from the histmakers: its OpenBLAS and the one ROOT pulls in end up in the
    same process and the least-squares solve inside the interpolator segfaults.
    The tridiagonal system for the knot slopes is small enough to solve with a
    plain Thomas sweep, which touches no BLAS at all.

    xnew is clamped to the range of x rather than extrapolated.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n = x.size - 1
    h = np.diff(x)
    slope = np.diff(y) / h

    # system for the slopes m_i = s'(x_i): C2 continuity in the interior rows,
    # continuous third derivative across the first and last interior knot
    lower, diag, upper, rhs = (np.zeros(n + 1) for _ in range(4))
    diag[1:n] = 2.0 * (h[:-1] + h[1:])
    lower[1:n] = h[1:]
    upper[1:n] = h[:-1]
    rhs[1:n] = 3.0 * (h[1:] * slope[:-1] + h[:-1] * slope[1:])

    dx0 = x[2] - x[0]
    diag[0] = h[1]
    upper[0] = dx0
    rhs[0] = ((h[0] + 2.0 * dx0) * h[1] * slope[0] + h[0] ** 2 * slope[1]) / dx0

    dxn = x[n] - x[n - 2]
    lower[n] = dxn
    diag[n] = h[n - 2]
    rhs[n] = (
        h[n - 1] ** 2 * slope[n - 2] + (2.0 * dxn + h[n - 1]) * h[n - 2] * slope[n - 1]
    ) / dxn

    for i in range(1, n + 1):
        w = lower[i] / diag[i - 1]
        diag[i] -= w * upper[i - 1]
        rhs[i] -= w * rhs[i - 1]
    m = np.empty(n + 1)
    m[n] = rhs[n] / diag[n]
    for i in range(n - 1, -1, -1):
        m[i] = (rhs[i] - upper[i] * m[i + 1]) / diag[i]

    xc = np.clip(xnew, x[0], x[-1])
    i = np.clip(np.searchsorted(x, xc, side="right") - 1, 0, n - 1)
    t = xc - x[i]
    c2 = (3.0 * slope[i] - 2.0 * m[i] - m[i + 1]) / h[i]
    c3 = (m[i] + m[i + 1] - 2.0 * slope[i]) / h[i] ** 2
    return y[i] + t * (m[i] + t * (c2 + t * c3))


def make_cvh_antiveto_helper(era=None):
    """Helper for the CVH correction to the *anti*-veto efficiency.

    Companion of define_cvh_weight: that one corrects the muons that have to be
    reconstructed, this one the muon that has to be missed. A second muon whose
    refit fails is absent from 'vetoMuons', so in data (only) a dimuon event can
    end up in the single-muon selection; the correction is an extra factor on
    the anti-veto SF applied to the unmatched postFSR gen muon.

    The factor is 1 + eps_data*(1 - SF_cvh)/(1 - eps_data), reconstructing
    1 - eps_data as antivetoSF*(1 - eps_MC) from the MC-truth veto efficiency
    map, which is the same input the anti-veto SF itself was derived from in
    scripts/analysisTools/w_mass_13TeV/smoothVetoSF.py. See the header.

    The map is stored in 2 GeV pt bins but 1 - eps_MC sits in the denominator
    and varies fast through the turn-on, so it is spline interpolated vs pt in
    each eta bin first, as smoothVetoSF.py does before building the SF. Using
    the raw bins instead would misestimate 1 - eps_MC by up to ~30% around
    20 GeV, hence the correction by as much.
    """
    eradir = era if era in ["2017", "2018"] else ""
    filename = (
        f"{data_dir}/muonSF/{eradir}/veto_global_SF/vetoEfficienciesEtaPt.pkl.lz4"
    )
    if not os.path.isfile(filename):
        # no such input for era 2017, where the veto SFs themselves are missing
        raise IOError(
            f"Couldn't read MC veto efficiency file {filename}, needed by the CVH"
            f" anti-veto correction for era {era}. Run with '--cvhBadModules veto'"
            " (drops the affected phase space instead) or 'none'."
        )
    logger.info(f"CVH anti-veto correction: MC veto efficiency read from {filename}")
    with lz4.frame.open(filename) as feff:
        alleff = pickle.load(feff)

    # categorical axes in python bindings always have an overflow bin, so use a
    # regular axis for the charge, as everywhere else in the SF helpers
    axis_charge = hist.axis.Regular(
        2, -2.0, 2.0, underflow=False, overflow=False, name="SF charge"
    )
    charges = {-1.0: "minus", 1.0: "plus"}

    heff = None
    for charge, charge_tag in charges.items():
        heff_charge = alleff[f"Wmunu_MC_eff_veto{charge_tag}_etapt"]
        axis_eta = heff_charge.axes[0]  # 48 bins in [-2.4, 2.4]
        axis_pt_in = heff_charge.axes[1]  # 27 bins of 2 GeV in [15, 69]
        if heff is None:
            axis_pt = hist.axis.Regular(
                5 * axis_pt_in.size,
                axis_pt_in.edges[0],
                axis_pt_in.edges[-1],
                underflow=False,
                overflow=False,
                name=axis_pt_in.name,
            )
            heff = hist.Hist(
                axis_eta,
                axis_pt,
                axis_charge,
                name="cvh_antiveto_effMC",
                storage=hist.storage.Weight(),
            )
        # nan would mean no event was selected in that bin, treat it as
        # efficiency 0 as smoothVetoSF.py does (none in the current input)
        values = np.nan_to_num(heff_charge.values(), nan=0.0)
        for ieta in range(axis_eta.size):
            heff.values()[ieta, :, axis_charge.index(charge)] = np.clip(
                _cubic_spline(axis_pt_in.centers, values[ieta], axis_pt.centers),
                0.0,
                1.0,
            )

    eff_pyroot = narf.hist_to_pyroot_boost(heff)
    helper = ROOT.wrem.cvh_antiveto_helper[type(eff_pyroot)](
        ROOT.std.move(eff_pyroot),
    )
    return helper


def define_cvh_antiveto_weight(
    df,
    helper,
    pt_col,
    eta_col,
    phi_col,
    charge_col,
    antiveto_sf_col=None,
    out_col="weight_cvhAntiVetoSF",
):
    """Define the extra anti-veto weight from make_cvh_antiveto_helper.

    The columns are those of the postFSR gen muon that failed the veto (charge
    = -99 when there is none, giving a weight of 1). antiveto_sf_col is the
    anti-veto SF column already folded into the event weight for that muon;
    pass None when no veto SF is applied, in which case the correction is taken
    with respect to the MC efficiency directly.

    Must be defined *after* the anti-veto SF it multiplies, and only on the MC
    samples that receive it.
    """
    if antiveto_sf_col is None:
        antiveto_sf_col = f"{out_col}_unitSF"
        df = df.Define(antiveto_sf_col, "1.0")
    df = df.Define(
        out_col, helper, [pt_col, eta_col, phi_col, charge_col, antiveto_sf_col]
    )
    return df, out_col
