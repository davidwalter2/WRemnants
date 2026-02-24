import hist
import numpy as np
from scipy import ndimage

from wremnants.utilities import common
from wums import boostHistHelpers as hh
from wums import logging

logger = logging.child_logger(__name__)


extended_pdf_datasets = [
    x for x in common.vprocs if not any(y in x for y in ["NNLOPS", "MiNLO"])
]


def replace_by_neighbors(vals, replace):
    if np.count_nonzero(replace) == vals.size:
        raise ValueError("Cannot replace all values with nearest non-zero neighbour")

    indices = ndimage.distance_transform_edt(
        replace, return_distances=False, return_indices=True
    )
    return vals[tuple(indices)]


def qcdByHelicityLabels():
    coeffs = ["const"] + [f"a{i}" for i in range(8)]
    scaleVars = ["muRmuF", "muR", "muF"]
    return [
        f"{var}_{coeff}{t}"
        for var in scaleVars
        for t in ["Up", "Down"]
        for coeff in coeffs
    ]


def qcdScaleNames():
    # Exclude central and extreme variations
    shifts = [
        "muRmuFDown",
        "muRDown",
        "",
        "muFDown",
        "",
        "muFUp",
        "muRUp",
        "",
        "muRmuFUp",
    ]
    return ["_".join(["QCDscale", s]) if s != "" else s for s in shifts]


def pdfNames(cardTool, pdf, skipFirst=True):
    size = 101
    names = cardTool.mirrorNames(f"pdf{{i}}{pdf}", size)
    if skipFirst:
        names[0] = ""
        names[size] = ""
    # TODO: This is probably not needed anymore, check with low PU
    if False and pdf == "NNPDF31":
        names[size - 2] = "pdfAlphas002Up"
        names[size - 1] = "pdfAlphas002Down"
        # Drop the mirrored alphaS variations
        names[size * 2 - 2] = ""
        names[size * 2 - 1] = ""
    return names


def pdfNamesAsymHessian(entries, pdfset=""):
    pdfNames = ["pdf0" + pdfset.replace("pdf", "")]
    if pdfset == "pdfHERAPDF20ext":
        entries -= 3
    pdfNames.extend(
        [
            f"pdf{int((j+2)/2)}{pdfset.replace('pdf', '')}{'Up' if j % 2 else 'Down'}"
            for j in range(entries - 1)
        ]
    )
    if pdfset == "pdfHERAPDF20ext":
        pdfNames.extend(
            [f"pdf{entries//2 + j}{pdfset.replace('pdf', '')}" for j in range(1, 4)]
        )
    return pdfNames


def pdfNamesSymHessian(entries, pdfset=""):
    return [f"pdf{i+1}{pdfset.replace('pdf', '')}" for i in range(entries)]


def pdfSymmetricShifts(hdiff, axis_name):
    sq = hh.multiplyHists(hdiff, hdiff)
    ss = sq[{axis_name: hist.sum}]
    rss = hh.sqrtHist(ss)
    return rss, rss


def pdfAsymmetricShifts(hdiff, axis_name):
    # Assuming that the last axis is the syst axis
    # TODO: add some check to verify this
    def shiftHist(vals, hdiff, axis_name):
        hnew = hdiff[{axis_name: 0}]
        vals = vals * vals
        hnew[...] = np.sum(vals, axis=-1)
        return hh.sqrtHist(hnew)

    ax = hdiff.axes[axis_name]
    underflow = hdiff.axes[axis_name].traits.underflow
    if type(ax) == hist.axis.StrCategory and all(
        ["Up" in x or "Down" in x for x in ax][1:]
    ):
        # Remove the overflow from the categorical axis
        end = int((ax.size - 1) / 2)
        upvals = hdiff[{axis_name: [x for x in ax if "Up" in x]}].values(flow=True)[
            ..., :end
        ]
        downvals = hdiff[{axis_name: [x for x in ax if "Down" in x]}].values(flow=True)[
            ..., :end
        ]
        if upvals.shape != downvals.shape:
            raise ValueError(
                "Malformed PDF uncertainty hist! Expect equal number of up and down vars"
            )
    else:
        end = ax.size + underflow
        upvals = hdiff.values(flow=True)[..., 1 + underflow : end : 2]
        downvals = hdiff.values(flow=True)[..., 2 + underflow : end : 2]

    # The error sets are ordered up,down,up,down...
    upshift = shiftHist(upvals, hdiff, axis_name)
    downshift = shiftHist(downvals, hdiff, axis_name)
    return upshift, downshift


def hessianPdfUnc(h, axis_name="pdfVar", uncType="symHessian", scale=1.0):
    symmetric = uncType == "symHessian"
    diff = hh.addHists(h, -1 * h[{axis_name: 0}]) * scale
    if diff.axes[axis_name].traits.overflow:
        diff[..., hist.overflow] = np.zeros_like(diff[{axis_name: 0}].view(flow=True))
    shiftFunc = pdfSymmetricShifts if symmetric else pdfAsymmetricShifts
    rssUp, rssDown = shiftFunc(diff, axis_name)
    hUp = hh.addHists(h[{axis_name: 0}], 1 * rssUp)
    hDown = hh.addHists(h[{axis_name: 0}], -1 * rssDown)
    return hUp, hDown
