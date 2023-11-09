import mplhep as hep
import matplotlib.pyplot as plt

from utilities import logging, common
from utilities.io_tools import output_tools
from utilities.io_tools.combinetf_input import get_fitresult
from wremnants import plot_tools

import pandas as pd
import numpy as np

import pdb

def get_mass_obs(filename):
    fitresult = get_fitresult(filename)

    val = fitresult["nois_outvals"][...][0]
    err = fitresult["nois_outcov"][...][0,0]**0.5

    return val, err

hep.style.use(hep.style.ROOT)

parser = common.plot_parser()
parser.add_argument("inputs", nargs="+", type=str, help="Paths to fitresult files")
args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)
outdir = output_tools.make_plot_dir(args.outpath, args.outfolder)

df = pd.DataFrame(args.inputs, columns=["path"])

df["base_name"] = df["path"].apply(lambda x: x.split("/")[-2])
df["analysis"] = df["base_name"].apply(lambda x: x.split("_")[0] )
df["mass"] = df["base_name"].apply(lambda x: int(x.split("_")[-1].split("MeV")[0].replace("massShift","")[1:]) * (-1 if x.endswith("Down") else 1) )
df[["mass_obs", "mass_err"]] = df["path"].apply(get_mass_obs).apply(pd.Series)

### make plot

pulls=True
diffs=False

for analysis, df_proc in df.groupby("analysis"):

    df_proc = df_proc.sort_values("mass")

    xarr = df_proc["mass"].values
    yarr = df_proc["mass_obs"].values * 100
    yerr = df_proc["mass_err"].values * 100

    ylim = min(yarr-yerr), max(yarr+yerr)

    if pulls:
        # pulls
        rlabel = "Pulls" #"(Meas.-True)/err(Meas.)"
        rarr = (xarr - yarr)/yerr
        rerr = np.zeros(len(yarr))
    elif diffs:
        # ratios
        rlabel = "Diff." #"(Meas.-True)/err(Meas.)"
        rarr = yarr - xarr
        rerr = yarr

    # extend ratio range by 20%
    rrange = min(rarr-rerr), max(rarr+rerr)
    rw = rrange[1] - rrange[0]
    rrange = rrange[0]-rw*0.1, rrange[1]+rw*0.1

    # extend x-axis range by 2%
    xlim=(min(xarr)-10, max(xarr)+10)
    xrange = xlim[1] - xlim[0]
    xlim = xlim[0]-xrange*0.01, xlim[1]+xrange*0.01

    fig, ax1, ax2 = plot_tools.figureWithRatio(xarr, "True mass", "Measured mass", ylim, rlabel, rrange, xlim=xlim, width_scale=1)

    ax1.plot([xlim[0], xlim[-1]], [xlim[0], xlim[-1]], linestyle="--", color="grey", label="Expectation")
    ax1.errorbar(xarr, yarr, yerr=yerr, marker=".", linestyle="", color="k", label="Measurement")

    ax2.plot([xlim[0], xlim[-1]], [0, 0], marker=".", linestyle="--", color="grey", label="Expectation")

    if pulls:
        ax2.plot(xarr, rarr,  linestyle="", marker=".", color="k")
    elif diffs:
        ax2.errorbar(xarr, rarr,  yerr=rerr, linestyle="", marker=".", color="k")

    plot_tools.addLegend(ax1, ncols=1, text_size=12,loc='upper left')

    scale = max(1, np.divide(*ax1.get_figure().get_size_inches()))
    hep.cms.label(ax=ax1, lumi=float(f"{args.lumi:.3g}"), fontsize=20*args.scaleleg*scale, 
        label=args.cmsDecor, data="True")

    outfile = f"massbias_{analysis}"
    if args.postfix:
        outfile += f"_{args.postfix}"
    plot_tools.save_pdf_and_png(outdir, outfile)
    plot_tools.write_index_and_log(outdir, outfile, args=args)
