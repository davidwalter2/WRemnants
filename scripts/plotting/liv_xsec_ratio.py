"""
Plot sigma_LIV / sigma_SM of the double differential SME model
(rabbit/param_models/sme_double_differential.py) for one Wilson coefficient,
either along mll for a fixed rapidity bin or along rapidity for a fixed mll bin.

The ratio is computed exactly as in LIV_rapidity.compute:
    sigma_LIV / sigma_SM = 1 + value * scale * sum_k cL_k(t) * resp_k / sigma_SM
where resp_k is the derivative of the bin integrated cross section w.r.t. the
lab frame coefficient k (c00, c33, c03) of the coefficient in the SU(2) basis of
arXiv:2011.02632 (c_u, c_d, d_u, c_c, c_s, d_c; LIV_rapidity SU2_BASIS: e.g. c_u shifts u with
c_uL = c_uR and d with c_dL = -c_dR), and cL_k(t) the lab frame coefficients
at sidereal hour t. The ratio oscillates with t, and the plotted value is its
maximum over a full sidereal day, taken separately in each bin.

With --vs mll and --smeModel, the same coefficient is also drawn as computed by
the mll-only model (rabbit/param_models/sme_model.py), from its precomputed
[time, mll] templates (integrated over the full rapidity range, 24 time bins):
    sigma_LIV / sigma_SM = 1 + value * scale * (L +- R) / 2 / sigma_SM
with + for c and - for d coefficients.

Examples:
    python scripts/plotting/liv_xsec_ratio.py cxxu --vs mll --fixedBin 3
    python scripts/plotting/liv_xsec_ratio.py dxzu --vs y --fixedRange 76 106 --value 5
    python scripts/plotting/liv_xsec_ratio.py cxxu --vs mll --inclusiveY --smeModel
    python scripts/plotting/liv_xsec_ratio.py cxxu --vs mll --fixedBin 3 --fullY --mllStep 0.5 -j 16
    python scripts/plotting/liv_xsec_ratio.py cxxu --vs mll --fixedBin 3 --smeModel --fullY --mllStep 0.5 -j 16
    python scripts/plotting/liv_xsec_ratio.py cxxu --vs y --fixedRange 76 106 --yStep 0.1 -j 16

With --mllStep, sme_model has no templates for the fine mll bins, and they are computed
here with the sme_model code (rabbit/param_models/liv/sme_functions.py), integrating
dsigma/dmll^2 over each fine bin with --smeModelPoints points (Simpson).

With --fullY, the double differential model is also drawn integrated over the full
kinematic rapidity range, |y| < ln(sqrt(S)/mll), which is what sme_model is integrated
over as well, so the two should agree.
"""

import argparse
import hashlib
import os
import pickle
from multiprocessing import Pool

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rabbit.param_models.liv_rapidity import functions as liv_functions
from rabbit.param_models.liv_rapidity.basis import (
    LIV_QUARKS,
    PREV_LIV_QUARKS,
    SM_QUARKS,
    bin_response,
    quark_response,
    sm_xsec,
)
from rabbit.param_models.liv_rapidity.functions import (
    S,
    lab_coeffs,
    make_matrix,
)
from rabbit.param_models.sme_double_differential import SU2_BASIS
from wums import output_tools, plot_tools

# same conventions as LIV_rapidity (sme_double_differential.py): xx is c_XX - c_YY, and
# xy, xz, yz have the signs of c^{XY}, c^{XZ}, c^{YZ} in Eq. 30 of arXiv:2011.02632
COEFFS = {
    "xx": 0.5 * (make_matrix(1, 1, 1.0) - make_matrix(2, 2, 1.0)),
    "xy": make_matrix(1, 2, 1.0),
    "xz": make_matrix(1, 3, 1.0),
    "yz": make_matrix(2, 3, 1.0),
    "ty": make_matrix(0, 2, -1.0),
    "tx": make_matrix(0, 1, -1.0),
}
# The coefficient is in the SU(2) basis of LIV_rapidity (SU2_BASIS): c_u, c_d, d_u, c_c, c_s, d_c (d_s an alias of d_c),
# each a fixed combination of the chiral quark coefficients c_qL, c_qR.
# sme_model (liv/sme_constants.py CL1-CL4) still uses the old conventions: xx with
# c_XX = -c_YY = 1 (i.e. c_XX - c_YY = 2) and the opposite sign for xy, xz, yz. Its curves
# are scaled by these factors, so that they stay comparable to the double differential at
# the same --value, as before the conventions of LIV_rapidity were changed.
SME_MODEL_CONV = {"xx": 0.5, "xy": -1.0, "xz": -1.0, "yz": -1.0}
Q_ID = {"d": 1, "u": 2, "s": 3, "c": 4, "b": 5}
MODE = {"c": 0, "d": 1}
EPS = 1e-6
PDF_NAME = "NNPDF31_nnlo_as_0118"

# default binning of scripts/histmakers/mz_dilepton_liv_rapidity.py
MLL_EDGES = [15, 30, 40, 45, 50, 55, 60, 65, 70, 76, 106, 110, 115, 120]
Y_EDGES = [-2.4, -1.4, -0.7, 0, 0.7, 1.4, 2.4]


def bin_cache_file(cache_dir, key, liv_quarks, with_sm_quarks=True):
    """per bin cache of the response basis of liv_quarks, with the SM cross section summed
    over SM_QUARKS (with_sm_quarks) or over liv_quarks only (oldest caches)"""
    tag_key = (key, list(liv_quarks), EPS, PDF_NAME) + (
        (SM_QUARKS,) if with_sm_quarks else ()
    )
    tag = hashlib.md5(repr(tag_key).encode()).hexdigest()[:12]
    return os.path.join(cache_dir, f"liv_bin_{tag}.npz")


def cached_bin_basis(args):
    """SM cross section summed over SM_QUARKS and response basis [mode, direction, quark]
    of one (mll, rapidity) bin, as in LIV_rapidity (basis.bin_response)"""
    key, cache_dir = args
    cache_file = bin_cache_file(cache_dir, key, LIV_QUARKS)
    if os.path.exists(cache_file):
        cached = np.load(cache_file)
        return float(cached["sm_total"]), cached["resp"]
    # reuse the responses of PREV_LIV_QUARKS (cached before charm was added)
    for with_sm_quarks in (True, False):
        old_file = bin_cache_file(cache_dir, key, PREV_LIV_QUARKS, with_sm_quarks)
        if not os.path.exists(old_file):
            continue
        cached = np.load(old_file)
        sm_total = float(cached["sm_total"])
        if not with_sm_quarks:
            sm_total += sm_xsec(key, [q for q in SM_QUARKS if q not in PREV_LIV_QUARKS])
        new_quarks = [q for q in LIV_QUARKS if q not in PREV_LIV_QUARKS]
        print(
            f"Integrating quarks {new_quarks} of bin mll=[{key[0]}, {key[1]}], y=[{key[2]}, {key[3]}]"
        )
        _, new_resp = quark_response(key, new_quarks)
        resp = np.zeros((2, 3, len(LIV_QUARKS)))
        for iq, q in enumerate(PREV_LIV_QUARKS):
            resp[:, :, LIV_QUARKS.index(q)] = cached["resp"][:, :, iq]
        for iq, q in enumerate(new_quarks):
            resp[:, :, LIV_QUARKS.index(q)] = new_resp[:, :, iq]
        break
    else:
        print(f"Integrating bin mll=[{key[0]}, {key[1]}], y=[{key[2]}, {key[3]}]")
        sm_total, resp = bin_response(key)
    np.savez(cache_file, sm_total=sm_total, resp=resp, key=np.array(key))
    return sm_total, resp


def sme_model_ratio(coeff, mll_edges, value, scale, liv_dir):
    """sigma_LIV / sigma_SM [time, mll] of the coefficient as in LIV (sme_model.py)"""
    if coeff[1:3] not in ["xx", "xy", "xz", "yz"]:
        raise ValueError(
            f"sme_model has no templates for '{coeff}', only c/d + xx, xy, xz, yz"
        )
    pdir = os.path.join(liv_dir, "precomputed_sigma")
    tag = f"{int(mll_edges[0])}_to_{int(mll_edges[-1])}_GeV_{len(mll_edges) - 1}_bins"
    with open(os.path.join(pdir, f"SM_{tag}.pkl"), "rb") as f:
        sm = pickle.load(f)
    if not np.allclose(sm["bins"], mll_edges):
        raise ValueError(
            f"sme_model templates have mll bins {sm['bins']}, plotting {list(mll_edges)}"
        )
    # sum of the chiral quark templates of the SU(2) basis coefficient, with the factor 1/2
    # of sme_model (1/2 (L +- R) for a single quark)
    shift = 0.0
    for (q, chirality), w in SU2_BASIS[(coeff[0], coeff[3])].items():
        with open(
            os.path.join(pdir, f"summation_{tag}_c{coeff[1:3]}_{q}_{chirality}.pkl"),
            "rb",
        ) as f:
            shift = shift + w * np.array(pickle.load(f)["values"])  # [time, mll]
    value = value * SME_MODEL_CONV[coeff[1:3]]
    return 1.0 + value * scale * 0.5 * shift / np.array(sm["values"])


def sme_model_bin(args):
    """SM cross section summed over quarks [] and the dsigma/dmll^2 integrals
    [quark (u, d, s), (pipi_L, pipj_L, pipi_R, pipj_R)] of one mll bin, in pb, with the
    sme_model code, as in liv/sme_precalculate.py"""
    (Qlo, Qhi), npts, cache_dir = args
    tag = hashlib.md5(
        repr(("sme_model", Qlo, Qhi, npts, PDF_NAME, "sm_quarks")).encode()
    ).hexdigest()[:12]
    cache_file = os.path.join(cache_dir, f"sme_model_bin_{tag}.npz")
    if os.path.exists(cache_file):
        cached = np.load(cache_file)
        return float(cached["sm_total"]), cached["integrals"]
    from scipy.integrate import simpson

    from rabbit.param_models.liv import sme_functions as sf

    print(f"Integrating sme_model bin mll=[{Qlo}, {Qhi}]")
    q2 = np.linspace(Qlo**2, Qhi**2, npts)
    vals = np.zeros(
        (npts, len(sf.quark_couplings), 5)
    )  # [Q2, quark, (sm, pipi, pipj, sum_L, sum_R)]
    sm_vals = np.zeros(
        npts
    )  # SM summed over all quarks (sm_quark_couplings), as SM_*.pkl
    for i, Q2 in enumerate(q2):
        sm_vals[i] = sf.d_sigma_sm(Q2, sf.sm_quark_couplings, precomputed=True)[0].sum()
        for iq, (flavor, e_f, g_fR, g_fL) in enumerate(sf.quark_couplings):
            vals[i, iq, 1:3] = sf.integrate_sigma_hat_prime_sme(Q2 / sf.s, flavor, Q2)
            vals[i, iq, 3] = sf.summation_terms(Q2, e_f, g_fL)
            vals[i, iq, 4] = sf.summation_terms(Q2, e_f, g_fR)
    GeV_to_pb = sf.factor * 0.389379 * 1e9
    sm_total = GeV_to_pb * simpson(sm_vals, x=q2)
    integrals = GeV_to_pb * np.stack(
        [
            simpson(vals[:, :, ip] * vals[:, :, ich], x=q2, axis=0)
            for ich in (3, 4)
            for ip in (1, 2)
        ],
        axis=-1,
    )
    np.savez(
        cache_file,
        sm_total=sm_total,
        integrals=integrals,
        key=np.array([Qlo, Qhi, npts]),
    )
    return sm_total, integrals


def sme_model_beams(hours):
    """contravariant beam momenta p1, p2 [hour, 4] in the SCF at sidereal hours, as in
    liv/sme_constants.py scf_beams (horizontal beam, rotation angle 2 pi h / 24 of LIV_rapidity)
    """
    from rabbit.param_models.liv.sme_constants import scf_beams

    pm, pn = scf_beams(hours)
    return np.array([p.numpy() for p in pm], dtype=float), np.array(
        [p.numpy() for p in pn], dtype=float
    )


def sme_model_fine_ratio(coeff, mll_edges, value, scale, hours, npts, cache_dir, jobs):
    """sigma_LIV / sigma_SM [hour, mll] of the coefficient as in LIV (sme_model.py), for
    mll bins without precomputed templates"""
    from rabbit.param_models.liv.sme_constants import CL1, CL2, CL3, CL4

    tensors = {"xx": CL1, "xy": CL2, "xz": CL3, "yz": CL4}
    if coeff[1:3] not in tensors:
        raise ValueError(
            f"sme_model has no templates for '{coeff}', only c/d + xx, xy, xz, yz"
        )
    tasks = [
        ((float(lo), float(hi)), npts, cache_dir)
        for lo, hi in zip(mll_edges[:-1], mll_edges[1:])
    ]
    if jobs > 1:
        with Pool(jobs) as pool:
            results = pool.map(sme_model_bin, tasks)
    else:
        results = [sme_model_bin(t) for t in tasks]
    sm_total = np.array([r[0] for r in results])  # [mll]
    C = np.array(tensors[coeff[1:3]].numpy(), dtype=float)
    pm, pn = sme_model_beams(hours)
    pipi = np.einsum("mn,tm,tn->t", C, pm, pm) + np.einsum("mn,tm,tn->t", C, pn, pn)
    pipj = np.einsum("mn,tm,tn->t", C, pm, pn) + np.einsum("mn,tm,tn->t", C, pn, pm)
    # sum of the chiral quark terms of the SU(2) basis coefficient
    shift = 0.0
    for (q, chirality), w in SU2_BASIS[(coeff[0], coeff[3])].items():
        iq = ["u", "d", "s"].index(q)  # order of sme_constants.quark_couplings
        integrals = np.array(
            [r[1][iq] for r in results]
        )  # [mll, (pipi_L, pipj_L, pipi_R, pipj_R)]
        i0 = 0 if chirality == "L" else 2
        shift = shift + w * (
            np.outer(pipi, integrals[:, i0]) + np.outer(pipj, integrals[:, i0 + 1])
        )
    value = value * SME_MODEL_CONV[coeff[1:3]]
    return 1.0 + value * scale * 0.5 * shift / sm_total


def sci_latex(x):
    """x as LaTeX scientific notation, e.g. 1e-06 -> 1 \\times 10^{-6}"""
    if x == 0:
        return "0"
    exp = int(np.floor(np.log10(abs(x))))
    mant = float(f"{x / 10**exp:.3g}")
    if abs(mant) >= 10:
        mant, exp = mant / 10, exp + 1
    if exp == 0:
        return f"{mant:g}"
    return rf"{mant:g} \times 10^{{{exp}}}"


def parse_coeff(coeff):
    """lab frame tensor and {(quark index, chirality): weight} of an SU(2) basis coefficient"""
    if (
        len(coeff) != 4
        or coeff[1:3] not in COEFFS
        or (coeff[0], coeff[3]) not in SU2_BASIS
    ):
        raise ValueError(
            f"Invalid coefficient '{coeff}', expected e.g. cxxu, dxzu "
            f"(SU(2) basis {[m + '..' + q for m, q in SU2_BASIS]}, direction {list(COEFFS)})"
        )
    weights = {
        (LIV_QUARKS.index(Q_ID[q]), chirality): w
        for (q, chirality), w in SU2_BASIS[(coeff[0], coeff[3])].items()
    }
    return COEFFS[coeff[1:3]], weights


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "coeff",
        help="Wilson coefficient in the SU(2) basis (c..u, c..d, d..u, c..s), e.g. cxxu, dxzu, ctyd",
    )
    parser.add_argument(
        "--vs", choices=["mll", "y"], default="mll", help="Variable on the x axis"
    )
    fixed = parser.add_mutually_exclusive_group(required=True)
    fixed.add_argument(
        "--fixedBin",
        type=int,
        help="Index of the bin of the other variable (rapidity bin if --vs mll, mll bin if --vs y)",
    )
    fixed.add_argument(
        "--fixedRange",
        type=float,
        nargs=2,
        metavar=("LO", "HI"),
        help="Range of the other variable",
    )
    fixed.add_argument(
        "--inclusiveY",
        action="store_true",
        help="Integrate over the full rapidity range of --yEdges (-2.4 < y < 2.4), only with --vs mll",
    )
    parser.add_argument(
        "--smeModel",
        action="store_true",
        help="Also draw the coefficient as computed by sme_model (only with --vs mll)",
    )
    parser.add_argument(
        "--fullY",
        action="store_true",
        help="Also draw the double differential integrated over the full kinematic rapidity range |y| < ln(sqrt(S)/mll) (only with --vs mll)",
    )
    parser.add_argument(
        "--logDeviation",
        action="store_true",
        help="Plot max_t sigma_LIV/sigma_SM - 1 on a log y axis",
    )
    parser.add_argument(
        "--value",
        type=float,
        default=1.0,
        help="Coefficient value, in units of --scale",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=1e-6,
        help="Scale of the coefficient (as in LIV_rapidity)",
    )
    parser.add_argument(
        "--nTimes",
        type=int,
        default=1440,
        help="Number of time points over the sidereal day used to find the maximum",
    )
    parser.add_argument(
        "--mllEdges", type=float, nargs="+", default=MLL_EDGES, help="mll bin edges"
    )
    parser.add_argument(
        "--mllStep",
        type=float,
        default=None,
        help="Use fine mll bins of this width (GeV) between the first and last --mllEdges, only with --vs mll",
    )
    parser.add_argument(
        "--smeModelPoints",
        type=int,
        default=21,
        help="Number of mll^2 points per fine bin to integrate sme_model with (Simpson), only with --mllStep",
    )
    parser.add_argument(
        "--yEdges", type=float, nargs="+", default=Y_EDGES, help="Rapidity bin edges"
    )
    parser.add_argument(
        "--yStep",
        type=float,
        default=None,
        help="Use fine rapidity bins of this width between the first and last --yEdges and outside of them up to the kinematic limit, only with --vs y",
    )
    parser.add_argument(
        "--jobs",
        "-j",
        type=int,
        default=1,
        help="Number of processes for the integration",
    )
    parser.add_argument(
        "--cacheDir",
        default=os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(liv_functions.__file__))),
            "precomputed_sigma",
            "bins",
        ),
        help="Directory to cache the per bin integrals",
    )
    parser.add_argument("-o", "--outpath", default="./", help="Base path for output")
    parser.add_argument(
        "-p", "--postfix", default=None, help="Postfix to append on output file name"
    )
    parser.add_argument(
        "--eoscp",
        action="store_true",
        help="Override use of xrdcp and use the mount instead",
    )
    args = parser.parse_args()

    mat, weights = parse_coeff(args.coeff)

    mll_edges, y_edges = np.array(args.mllEdges), np.array(args.yEdges)
    if args.mllStep is not None:
        if args.vs != "mll":
            raise ValueError("--mllStep requires --vs mll")
        nbins = int(round((mll_edges[-1] - mll_edges[0]) / args.mllStep))
        mll_edges = np.linspace(mll_edges[0], mll_edges[-1], nbins + 1)
    if args.yStep is not None:
        if args.vs != "y":
            raise ValueError("--yStep requires --vs y")
        nbins = int(round((y_edges[-1] - y_edges[0]) / args.yStep))
        y_edges = np.linspace(y_edges[0], y_edges[-1], nbins + 1)
    scan_edges, other_edges = (
        (mll_edges, y_edges) if args.vs == "mll" else (y_edges, mll_edges)
    )
    if args.vs == "y" and (args.inclusiveY or args.smeModel or args.fullY):
        raise ValueError("--inclusiveY, --smeModel and --fullY require --vs mll")
    if args.inclusiveY:
        flo, fhi = y_edges[0], y_edges[-1]
    elif args.fixedRange is not None:
        flo, fhi = args.fixedRange
    else:
        if not 0 <= args.fixedBin < len(other_edges) - 1:
            raise ValueError(
                f"--fixedBin {args.fixedBin} out of range, {len(other_edges) - 1} bins: {other_edges}"
            )
        flo, fhi = other_edges[args.fixedBin], other_edges[args.fixedBin + 1]

    if args.vs == "y":
        # extend the rapidity range outside of the acceptance of --yEdges up to the kinematic
        # limit |y| < ln(sqrt(S)/mll) at the lower mll edge, in steps of --yStep (one bin
        # on each side without)
        y_kin = float(np.log(np.sqrt(S) / flo))
        outer = []
        for acc, sign in [(y_edges[-1], 1), (-y_edges[0], -1)]:
            if y_kin <= acc:
                outer.append(np.array([]))
                continue
            n = (
                1
                if args.yStep is None
                else max(1, int(round((y_kin - acc) / args.yStep)))
            )
            outer.append(sign * np.linspace(acc, y_kin, n + 1)[1:])
        scan_edges = np.concatenate([outer[1][::-1], y_edges, outer[0]])

    os.makedirs(args.cacheDir, exist_ok=True)
    hours = np.linspace(0, 24, args.nTimes, endpoint=False)
    cL = np.array(
        [[float(v) for v in lab_coeffs(h, mat)] for h in hours]
    )  # [hour, direction]

    def max_ratio(keys):
        """max over the sidereal day of sigma_LIV / sigma_SM in each (mll lo, mll hi, y lo, y hi) bin"""
        tasks = [(k, args.cacheDir) for k in keys]
        if args.jobs > 1:
            with Pool(args.jobs) as pool:
                results = pool.map(cached_bin_basis, tasks)
        else:
            results = [cached_bin_basis(t) for t in tasks]
        sm_total = np.array([r[0] for r in results])  # [bin]
        # SU(2) basis coefficient: sum of the chiral quark responses, as in LIV_rapidity
        # (r[1][mode] is the response to c_qL = +-c_qR = 1, mode c: +, d: -)
        resp = np.array(
            [
                sum(
                    w
                    * 0.5
                    * (r[1][0, :, iq] + (1 if ch == "L" else -1) * r[1][1, :, iq])
                    for (iq, ch), w in weights.items()
                )
                for r in results
            ]
        ).T  # [direction, bin]
        return (1.0 + args.value * args.scale * (cL @ resp) / sm_total).max(
            axis=0
        )  # [bin]

    keys = []
    for lo, hi in zip(scan_edges[:-1], scan_edges[1:]):
        mll, y = ((lo, hi), (flo, fhi)) if args.vs == "mll" else ((flo, fhi), (lo, hi))
        keys.append(tuple(float(v) for v in (*mll, *y)))

    bin_labels = [
        f"[{lo:g}, {hi:g}]" for lo, hi in zip(scan_edges[:-1], scan_edges[1:])
    ]
    # max_t sigma_LIV/sigma_SM - 1 per bin, also written to the log file
    max_deviation = {}

    # maximum of the sidereal modulation in each bin
    rmax = max_ratio(keys)
    max_deviation["sme_double_differential"] = {
        b: float(r - 1) for b, r in zip(bin_labels, rmax)
    }
    for b, r in zip(bin_labels, rmax):
        print(f"{b}: max sigma_LIV/sigma_SM - 1 = {r - 1:.4e}")

    fig, ax = plt.subplots(figsize=(8, 8))

    def draw(r, **kwargs):
        if args.logDeviation:
            # deviation from 1, non-positive values can not be shown on the log axis
            if np.any(r - 1 <= 0):
                print(
                    f"Warning: {kwargs['label']} has max_t sigma_LIV/sigma_SM - 1 <= 0 in some bins, not shown on the log axis"
                )
            r = np.where(r - 1 > 0, r - 1, np.nan)
        ax.stairs(r, scan_edges, baseline=None, linewidth=2, **kwargs)

    def draw_acceptance(r, inside_kwargs, outside_kwargs):
        """draw the bins inside and outside the rapidity acceptance of --yEdges separately"""
        centers = 0.5 * (scan_edges[:-1] + scan_edges[1:])
        inside = (centers > y_edges[0]) & (centers < y_edges[-1])
        if args.logDeviation:
            if np.any(r - 1 <= 0):
                print(
                    "Warning: max_t sigma_LIV/sigma_SM - 1 <= 0 in some bins, not shown on the log axis"
                )
            r = np.where(r - 1 > 0, r - 1, np.nan)
        # one stairs per contiguous run of bins inside or outside of the acceptance
        runs, start = [], 0
        for i in range(1, len(r) + 1):
            if i == len(r) or inside[i] != inside[start]:
                runs.append((start, i, bool(inside[start])))
                start = i
        # inside first, for the legend order
        labeled = set()
        for start, stop, is_inside in sorted(runs, key=lambda run: not run[2]):
            kwargs = dict(inside_kwargs if is_inside else outside_kwargs)
            if is_inside in labeled:
                kwargs.pop("label", None)
            labeled.add(is_inside)
            ax.stairs(
                r[start:stop],
                scan_edges[start : stop + 1],
                baseline=None,
                linewidth=2,
                **kwargs,
            )

    fixed_label = (
        rf"${flo:g} < y < {fhi:g}$"
        if args.vs == "mll"
        else rf"${flo:g} < m_{{\ell\ell}} < {fhi:g}$ GeV"
    )
    # same color for both curves of the coefficient, the rapidity binned one dashed
    color = "C0"
    if args.vs == "mll":
        draw(
            rmax,
            color=color,
            linestyle="--",
            label=rf"$d^2\sigma/(d(m_{{\ell\ell}})^2\,dy)$" + "\n" + fixed_label,
        )
    else:
        acc = rf"{y_edges[-1]:g}" if np.isclose(-y_edges[0], y_edges[-1]) else None
        draw_acceptance(
            rmax,
            dict(color=color, label=rf"$|y| < {acc}$" if acc else "inside acceptance"),
            dict(color="C1", label=rf"$|y| > {acc}$" if acc else "outside acceptance"),
        )
    if args.fullY:
        # double differential integrated over the full kinematic rapidity range, as sme_model
        ymax = [float(np.log(np.sqrt(S) / lo)) for lo in scan_edges[:-1]]
        rmax_full = max_ratio(
            [
                (float(lo), float(hi), -yk, yk)
                for lo, hi, yk in zip(scan_edges[:-1], scan_edges[1:], ymax)
            ]
        )
        max_deviation["sme_double_differential_fullY"] = {
            b: float(r - 1) for b, r in zip(bin_labels, rmax_full)
        }
        for b, r in zip(bin_labels, rmax_full):
            print(f"{b}: full y max sigma_LIV/sigma_SM - 1 = {r - 1:.4e}")
        # solid in the color of the coefficient, unless sme_model is drawn as well
        draw(
            rmax_full,
            color="C1" if args.smeModel else color,
            linestyle=":" if args.smeModel else "-",
            label=r"$d\sigma/d(m_{\ell\ell})^2$",
        )
    if args.smeModel:
        if args.mllStep is None:
            liv_dir = os.path.join(
                os.path.dirname(
                    os.path.dirname(os.path.abspath(liv_functions.__file__))
                ),
                "liv",
            )
            rmax_old = sme_model_ratio(
                args.coeff, scan_edges, args.value, args.scale, liv_dir
            ).max(axis=0)
        else:
            # no sme_model templates for the fine bins, compute them with the sme_model code
            rmax_old = sme_model_fine_ratio(
                args.coeff,
                scan_edges,
                args.value,
                args.scale,
                hours,
                args.smeModelPoints,
                args.cacheDir,
                args.jobs,
            ).max(axis=0)
        max_deviation["sme_model"] = {
            b: float(r - 1) for b, r in zip(bin_labels, rmax_old)
        }
        for b, r in zip(bin_labels, rmax_old):
            print(f"{b}: sme_model max sigma_LIV/sigma_SM - 1 = {r - 1:.4e}")
        # distinguish from the integrated double differential, which has the same label
        draw(
            rmax_old,
            color=color,
            label=r"$d\sigma/d(m_{\ell\ell})^2$"
            + (" (sme_model)" if args.fullY else ""),
        )
    if args.logDeviation:
        ax.set_yscale("log")
    else:
        ax.axhline(1.0, color="black", linestyle="--", linewidth=1)
    xlabel = r"$m_{\ell\ell}$ [GeV]" if args.vs == "mll" else r"$y$"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(
        r"$\sigma_{\mathrm{LIV}} / \sigma_{\mathrm{SM}}$"
        + (" $-$ 1" if args.logDeviation else "")
    )
    ax.set_xlim(scan_edges[0], scan_edges[-1])
    title = rf"${args.coeff[0]}^{{{args.coeff[1:3]}}}_{{{args.coeff[3]}}} = {sci_latex(args.value * args.scale)}$"
    if args.vs == "y":
        title += ", " + fixed_label
    ax.set_title(title)
    if not args.logDeviation:
        ax.ticklabel_format(axis="y", useOffset=False)
    ax.legend(loc="best")
    fig.tight_layout()

    outdir = output_tools.make_plot_dir(args.outpath, eoscp=args.eoscp)
    fixed_tag = f"{'y' if args.vs == 'mll' else 'mll'}_{flo:g}_{fhi:g}".replace(
        "-", "m"
    ).replace(".", "p")
    to_join = [
        "liv_ratio",
        args.coeff,
        f"vs_{args.vs}",
        fixed_tag,
        "smeModel" if args.smeModel else None,
        "fullY" if args.fullY else None,
        (
            f"mllStep_{args.mllStep:g}".replace(".", "p")
            if args.mllStep is not None
            else None
        ),
        f"yStep_{args.yStep:g}".replace(".", "p") if args.yStep is not None else None,
        "logDeviation" if args.logDeviation else None,
        args.postfix,
    ]
    outfile = "_".join(filter(lambda x: x, to_join))
    plot_tools.save_pdf_and_png(outdir, outfile, fig)
    output_tools.write_index_and_log(
        outdir,
        outfile,
        analysis_meta_info={"max_t sigma_LIV/sigma_SM - 1": max_deviation},
        args=args,
    )
    if output_tools.is_eosuser_path(args.outpath) and args.eoscp:
        output_tools.copy_to_eos(outdir, args.outpath)


if __name__ == "__main__":
    main()
