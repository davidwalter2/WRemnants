import argparse
import pickle

import hist
import lz4.frame
import numpy as np

from wremnants.postprocessing.datagroups.datagroups import Datagroups
from wums import logging

# Produce the MC efficiency (effMC) input for the in-situ muon efficiency
# method from the 4-category probe spectra emitted by mz_dilepton.py
# (effMCprobe_{nominal,failIso,failHLT,failID}, MC only).
#
#   nominal : 2HLT, probe passes ID & HLT & Iso
#   failIso : 2HLT, probe passes ID & HLT, fails Iso
#   failHLT : 1HLT, probe passes ID, fails HLT
#   failID  : 1HLT, probe fails ID
#
# Per-step MC efficiency (raw per-bin), measured from the Zmumu signal only:
#   eff_iso     =  nominal                       / (nominal + failIso)
#   eff_trigger = (nominal + failIso)            / (nominal + failIso + failHLT)
#   eff_idip    = (nominal + failIso + failHLT)  / (nominal + failIso + failHLT
#                                                              + failID)

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputFile", type=str, required=True)
parser.add_argument("-o", "--outputFile", type=str, default="insitu_effMC.pkl.lz4")
parser.add_argument(
    "--process",
    type=str,
    default="Zmumu",
    help="Datagroup to measure effMC from (in-situ uses the signal only)",
)
parser.add_argument("--debug", action="store_true", help="Print debug output")
args = parser.parse_args()

logger = logging.setup_logger("make_insitu_effMC", 4 if args.debug else 3)

categories = ["nominal", "failIso", "failHLT", "failID"]
steps = ["idip", "trigger", "iso"]

datagroups = Datagroups(args.inputFile)
if datagroups.mode != "z_dilepton":
    raise ValueError("Expected input is the output from the dilepton histmaker")

for cat in categories:
    datagroups.loadHistsForDatagroups(f"effMCprobe_{cat}", syst="")

groups = datagroups.groups
if args.process not in groups:
    raise ValueError(
        f"Process group '{args.process}' not found in input "
        f"(available: {list(groups.keys())})"
    )

# probe spectra for the requested signal group, summed weights per bin
h = {cat: groups[args.process].hists[f"effMCprobe_{cat}"] for cat in categories}
n = {cat: h[cat].values(flow=True) for cat in categories}

# template axes [eta, pt, charge] taken from the input hists
ref_axes = list(h["nominal"].axes)
axis_step = hist.axis.StrCategory(steps, name="insitu_step")
effMC = hist.Hist(*ref_axes, axis_step, name="effMC", storage=hist.storage.Double())


def _ratio(num, den):
    out = np.zeros_like(num, dtype=float)
    np.divide(num, den, out=out, where=den > 0)
    return out


pass2HLT = n["nominal"] + n["failIso"]
passID = pass2HLT + n["failHLT"]
allProbes = passID + n["failID"]

eff = {
    "iso": _ratio(n["nominal"], pass2HLT),
    "trigger": _ratio(pass2HLT, passID),
    "idip": _ratio(passID, allProbes),
}

for step in steps:
    effMC.view(flow=True)[..., axis_step.index(step)] = eff[step]
    e = eff[step]
    inrange = (e > 0) & (e < 1)
    nfilled = int(inrange.sum())
    if nfilled:
        sel = e[inrange]
        logger.info(
            f"effMC[{step}]: {nfilled} bins in (0,1), "
            f"mean={sel.mean():.4f} min={sel.min():.4f} max={sel.max():.4f}; "
            f"{int((e >= 1).sum())} bins ==1 (low stats), "
            f"{int((e <= 0).sum())} empty"
        )
    else:
        logger.warning(
            f"effMC[{step}]: no bins in (0,1) "
            f"({int((e >= 1).sum())} bins ==1, {int((e <= 0).sum())} empty) "
            f"- expected with very low statistics"
        )

with lz4.frame.open(args.outputFile, "wb") as fout:
    pickle.dump(
        {"effMC": effMC, "process": args.process, "steps": steps},
        fout,
        protocol=pickle.HIGHEST_PROTOCOL,
    )
logger.info(f"Wrote in-situ effMC to {args.outputFile}")
