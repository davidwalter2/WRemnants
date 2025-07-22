import matplotlib.pyplot as plt
from matplotlib import ticker

from utilities import parsing
from utilities.io_tools import hepdata_tools, rabbit_input
from wums import output_tools, plot_tools

parser = parsing.plot_parser()
parser.add_argument(
    "-r",
    "--reffile",
    required=True,
    type=str,
    help="Combine fitresult file for nominal result",
)
parser.add_argument("--print", action="store_true", help="Print results")
parser.add_argument(
    "--diffToCentral", action="store_true", help="Show difference to central result"
)
parser.add_argument(
    "--saveForHepdata",
    action="store_true",
    help="Save output as ROOT to prepare HEPData",
)
args = parser.parse_args()

isW = "WMass" in args.reffile
isWm = isW and "_Wm" in args.reffile
basename = args.reffile

postfix_names = [
    "",
    "_scetlib_dyturboN3p1LL",
    "_scetlib_nnlojetN3p1LLN3LO",
    "_scetlib_dyturboN4p0LL",  # "_dyturboN3LLp",
    "_scetlib_nnlojetN4p0LLN3LO",
    "_dataPtllRwgt",
]
labels = [
    "N$^{3{+}0}$LL+NNLO",
    "N$^{3{+}1}$LL+NNLO",
    "N$^{3{+}1}LL+N^{{3}}LO$",
    "N$^{4{+}0}$LL+NNLO",
    "N$^{4{+}0}LL+N^{{3}}LO$",
    r"$\mathit{p}_{T}^{\ell\ell}$ rwgt., N$^{3{+}0}$LL unc.",
]

if isW:
    postfix_names.append("_CombinedPtll")
    labels.append(
        r"Combined $\mathit{p}_{T}^{\ell\ell}$ fit," + "\n N$^{3{+}0}$LL unc."
    )

dfs = rabbit_input.read_all_groupunc_df(
    [args.reffile.format(postfix=p) for p in postfix_names],
    names=labels,
    uncs=["pTModeling"],
)

if isW:
    xlim = [80275, 80360] if isWm else [80331, 80372]
else:
    xlim = [91160, 91280] if "flipEvenOdd" not in basename else [91170, 91290]

if args.print:
    for k, v in dfs.iterrows():
        print(v.iloc[0], round(v.iloc[1], 1), round(v.iloc[3], 1), round(v.iloc[2], 2))

central = dfs.iloc[0, :]

xlabel = "".join(
    [r"$\mathit{m}_{", "W" if isW else "Z", "^{{-}}" if isWm else "", "}$ (MeV)"]
)

central_val = central["value"]
if args.diffToCentral:
    if args.saveForHepdata:
        # save also the original absolute value
        dfs["absolute_value"] = dfs["value"].values
    dfs["value"] -= central_val
    xlim = [xlim[0] - central_val, xlim[1] - central_val]
    central_val = 0
    xlabel = r"$\Delta$" + xlabel

fig = plot_tools.make_summary_plot(
    central_val,
    central["err_total"],
    central["err_pTModeling"],
    "N$^{3{+}0}$LL+NNLO\n (nominal)",
    dfs.iloc[1:, :],
    colors="auto",
    xlim=xlim,
    xlabel=xlabel,
    legend_loc="upper left",
    legtext_size="small",
    logoPos=0,
    cms_label=args.cmsDecor,
    lumi=16.8,
    padding=5,
)
ax = plt.gca()
minor_tick = 10 if isWm else 5
ax.yaxis.set_major_locator(ticker.MultipleLocator(2 * minor_tick))
ax.xaxis.set_major_locator(ticker.MultipleLocator(2 * minor_tick))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(minor_tick))
ax.xaxis.grid(False, which="both")
ax.yaxis.grid(False, which="both")

eoscp = output_tools.is_eosuser_path(args.outpath)
outdir = output_tools.make_plot_dir(args.outpath, args.outfolder, eoscp=eoscp)

outname = f"{'Wmass' if isW else 'Wlike'}_modeling_summary"
if args.postfix:
    outname += f"_{args.postfix}"

plot_tools.save_pdf_and_png(outdir, outname, fig)
output_tools.write_index_and_log(outdir, outname)

if args.saveForHepdata:
    column_labels = [xlabel, "Total uncertainty", "Model uncertainty"]
    if args.diffToCentral:
        column_labels.append(xlabel.replace(r"$\Delta$", ""))

    hepdata_tools.make_mass_summary_histogram(
        dfs, f"{outdir}/{outname}.root", column_labels
    )

if eoscp:
    output_tools.copy_to_eos(outdir, args.outpath, args.outfolder)
