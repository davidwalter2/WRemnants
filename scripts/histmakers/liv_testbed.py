import os

import hist

import narf
from wremnants.production import muon_efficiencies_liv
from wremnants.production.datasets.dataset_tools import getDatasets
from wremnants.production.histmaker_tools import (
    aggregate_groups,
    scale_to_data,
    write_analysis_output,
)
from wremnants.utilities import common, parsing
from wums import logging

analysis_label = common.analysis_label(os.path.basename(__file__))
parser, initargs = parsing.common_parser(analysis_label)

axis_eta = hist.axis.Variable(
    [-2.4, -1.40655, -0.68156, -0.00848, 0.66796, 1.4006, 2.4], name="eta_probe"
)


axis_pt_high = hist.axis.Variable(
    [
        15,
        # 21,
        25,
        32.35393,
        35.70991,
        38.30856,
        40.43642,
        42.22635,
        43.92092,
        45.87573,
        48.56281,
        53.1789,
        80,
    ],
    name="pt_tag",
)

axis_mll = hist.axis.Variable(
    [15, 30, 40, 45, 50, 55, 60, 65, 70, 76, 106, 110, 115, 120], name="mll"
)


redo_cdf = True
args = parser.parse_args()
logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)
era = args.era
calib_filepaths = common.calib_filepaths

muon_efficiency_helper = muon_efficiencies_liv.make_muon_efficiency_helper(era="2016BG")

print(muon_efficiency_helper({1.0, 2.0}, {30.0, 35.0}))

datasets = getDatasets(
    maxFiles=args.maxFiles,
    filt=args.filterProcs,
    excl=args.excludeProcs,
    nanoVersion="v9",
    base_path=args.dataPath,
    extended="msht20an3lo" not in args.pdfs,
    era=era,
)

########################################################


def build_graph(df, dataset):

    logger.info(f"fomrbuild graph for dataset: {dataset.name}")
    era = args.era
    results = []

    if not dataset.is_data:

        columnsForSF = [
            "Muon_eta",
            "Muon_pt",
        ]

        df = df.Define(
            "reco_scalefactor_weight",
            muon_efficiency_helper,
            columnsForSF,
        )

        df = df.Define("weight_expr", "genWeight * reco_scalefactor_weight")

    else:
        pass

    return results


logger.debug(f"Datasets are {[d.name for d in datasets]}")
resultdict = narf.build_and_run(datasets[::-1], build_graph)

if not args.noScaleToData:
    scale_to_data(resultdict)
    aggregate_groups(datasets, resultdict, args.aggregateGroups)

write_analysis_output(
    resultdict, f"{os.path.basename(__file__).replace('py', 'hdf5')}", args
)
