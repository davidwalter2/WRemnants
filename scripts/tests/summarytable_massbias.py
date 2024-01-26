import ROOT
import pandas as pd
import numpy as np
from scipy.stats import chi2

from utilities import logging, common
from utilities.io_tools import output_tools, tex_tools
from utilities.io_tools.combinetf_input import get_fitresult

import pdb

translate = {
    "WtaunuAsSig": r"$\mathrm{W}\rightarrow \tau\nu$ as sig.",
    "WtaunuAsBkg": r"$\mathrm{W}\rightarrow \tau\nu$ as bkg.",
    "nominal": r"-",
    "mtcut": r"$\Delta m^\mathrm{W}$",
    "dphiCut": r"$\Delta \phi$",
    "dphiCut_mtcut": r"$\Delta \phi\ \mathrm{and}\ \Delta m^\mathrm{W}$",
}


def read_fitresult(filename):
    try:
        rfile = ROOT.TFile.Open(filename)
        ttree = rfile.Get("fitresults")
        ttree.GetEntry(0)

        if hasattr(ttree,"massShiftW100MeV"):
            m = ttree.massShiftW100MeV*100
            merr = ttree.massShiftW100MeV_err*100
        else:
            m = 0
            merr = 0

        status = ttree.status
        errstatus = ttree.errstatus
        edmval = ttree.edmval

    except IOError as e:
        return 0, 1, 0, 0, 0, 0, 0, 0
        
    return status, errstatus, edmval, m, merr

parser = common.plot_parser()
parser.add_argument("inputs", nargs="+", type=str, help="Paths to fitresult files")
args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)
outdir = output_tools.make_plot_dir(args.outpath, args.outfolder)

df = pd.DataFrame(args.inputs, columns=["path"])

df[["status", "errstatus", "edmval", "mass_obs", "mass_err"]] = df["path"].apply(read_fitresult).apply(pd.Series)

df["name_parts"] = df["path"].apply(lambda x: [y for y in filter(lambda z: z, x.split("/"))])

df["dataset"] = df["name_parts"].apply(lambda x: x[-4].split("_")[0])
df["column_name"] = df["name_parts"].apply(lambda x: x[-3])

df["dataset"] = df["dataset"].apply(lambda x: translate.get(x, x))
df["column_name"] = df["column_name"].apply(lambda x: translate.get(x, x))

tex_tools.make_latex_table(df, output_dir=outdir, output_name=f"table_massVariation10", 
    column_title="Gen cuts", caption="Caption", 
    label="Pseudodata", sublabel="",
    column_name="column_name", row_name="dataset", 
    cell_columns=["mass_obs", "mass_err"], color_condition=lambda x, y: x-100 > 10., cell_format=lambda x, y: f"${round(x,2)} \pm {round(y,2)}"+" \mathrm{MeV}$")

# for channel, df_c in df.groupby("channel"):

    # tex_tools.make_latex_table(df_c, output_dir=outdir, output_name=f"table_status_{channel}", 
    #     column_title="Axes", caption="Fit status and error status.", label="Pseudodata", sublabel="",
    #     column_name="column_name", row_name="dataset", 
    #     cell_columns=["status", "errstatus"], color_condition=lambda x, y: y !=0, cell_format=lambda x, y: f"{int(x)} ({int(y)})")

    # tex_tools.make_latex_table(df_c, output_dir=outdir, output_name=f"table_edmval_{channel}", 
    #     column_title="Axes", caption="Estimated distance to minimum.", label="Pseudodata", sublabel="",
    #     column_name="column_name", row_name="dataset", 
    #     cell_columns=["edmval"], color_condition=lambda x: False, cell_format=lambda x: x)

    # tex_tools.make_latex_table(df_c, output_dir=outdir, output_name=f"table_mass_{channel}", 
    #     column_title="Axes", caption="Mass and uncertainty.", label="Pseudodata", sublabel="",
    #     column_name="column_name", row_name="dataset", 
    #     cell_columns=["mass_obs", "mass_err"], color_condition=lambda x, y: x > y, cell_format=lambda x, y: f"${round(x*100,2)}\, \pm {round(y*100,2)}$")

