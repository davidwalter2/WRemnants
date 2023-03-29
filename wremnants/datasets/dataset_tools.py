import narf
from utilities import logging
import subprocess
import glob

logger = logging.child_logger(__name__)

#set the debug level for logging incase of full printout 
from wremnants.datasets import files_v9
from wremnants.datasets import dataDictV8
from wremnants.datasets.datasetDict_gen import genDataDict
from wremnants.datasets.crosssections import xsec_13TeV

from wremnants.datasets.dataset_tools import filterProcs, excludeProcs, makeFilelist

logger = logging.child_logger(__name__)

lumicsv = f"{pathlib.Path(__file__).parent.parent}/data/bylsoutput.csv"
lumijson = f"{pathlib.Path(__file__).parent.parent}/data/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"

def getDatasets(maxFiles=-1, filt=None, excl=None, mode=None, base_path=None, nanoVersion="v9", 
        data_tag="TrackFitV722_NanoProdv2", mc_tag="TrackFitV718_NanoProdv1"):
    if not base_path:
        hostname = socket.gethostname()
        if hostname == "lxplus8s10.cern.ch":
            base_path = "/scratch/shared/NanoAOD"
        if hostname == "cmswmass2.cern.ch":
            base_path = "/data/shared/NanoAOD"
        elif "mit.edu" in hostname:
            base_path = "/scratch/submit/cms/wmass/NanoAOD"
        elif hostname == "cmsanalysis.pi.infn.it":
            base_path = "/scratchnvme/wmass/NANOV9/postVFP" #temporary

    logger.info(f"Loading samples from {base_path}.")

    if nanoVersion == "v9":
        dataDict = files_v9.files
    else:
        raise ValueError("Only NanoAODv8 and NanoAODv9 are supported")

    if mode == "gen":
        dataDict.update(genDataDict)

    narf_datasets = []
    for era, eraDict in dataDict.items():
        for sample, info in eraDict.items():
            if sample in genDataDict:
                base_path = base_path.replace("NanoAOD", "NanoGen")

            is_data = "data" in ["SingleMuon", ]

            prod_tag = data_tag if is_data else mc_tag 
            paths = makeFilelist(info["filepaths"], maxFiles, format_args=dict(BASE_PATH=base_path, NANO_PROD_TAG=prod_tag))

            if not paths:
                logger.warning(f"Failed to find any files for dataset {sample}. Looking at {info['filepaths']}. Skipping!")
                continue

            narf_info = dict(
                name=sample,
                filepaths=paths,
            )

            if is_data:
                if mode == "gen":
                    continue
                narf_info.update(dict(
                    is_data=True,
                    lumi_csv=lumicsv,
                    lumi_json=lumijson,
                ))
            else:
                xsec = xsec_13TeV.get(sample, None)
                if xsec is None:
                    logger.warning(f"No cross section found for sample {sample}!")

                narf_info.update(dict(
                    xsec=xsec,
                    )
                )
            narf_datasets.append(narf.Dataset(**narf_info))

    narf_datasets = filterSamples(filt, narf_datasets)
    narf_datasets = excludeSamples(excl, narf_datasets)
    
    for sample in narf_datasets:
        if not sample.filepaths:
            logger.warning(f"Failed to find any files for sample {sample.name}!")

    return narf_datasets



def buildXrdFileList(path, xrd):
    xrdpath = path[path.find('/store'):]
    logger.debug(f"Looking for path {xrdpath}")
    # xrdfs doesn't like wildcards, just use the mount if they are included
    if "*" not in path:
        f = subprocess.check_output(['xrdfs', f'root://{xrd}', 'ls', xrdpath]).decode(sys.stdout.encoding)
        return filter(lambda x: "root" in x[-4:], f.split())
    else:
        return [f"root://{xrd}/{f}" for f in glob.glob(path)]

#TODO add the rest of the samples!
def makeFilelist(paths, maxFiles=-1, format_args={}):
    filelist = []
    for path in paths:
        if format_args:
            path = path.format(**format_args)
            logger.debug(f"Reading files from path {path}")
        filelist.extend(glob.glob(path) if path[:4] != "/eos" else buildXrdFileList(path, "eoscms.cern.ch"))
    return filelist if maxFiles < 0 else filelist[:maxFiles]


def selectSample(selection, datasets):
    return list(filter(lambda x, s=selection: s in x.name, datasets))

def selectSamples(selections, datasets):
    new_datasets = []
    for selection in selections:
        new_datasets += selectSample(selection, datasets)

    # remove duplicates selected by multiple filters
    new_datasets = list(set(new_datasets))
    return new_datasets

def filterSamples(filters, datasets):
    if filters:
        if isinstance(filters, list):
            new_datasets = selectSamples(filters, datasets)
        elif isinstance(filters, str):
            new_datasets = selectSample(filters, datasets)
        else:
            new_datasets = list(filter(filters, datasets))
    else:
        return datasets

    if len(new_datasets) == 0:
        logger.warning("Try to filter samples but didn't find any match. Continue without filtering.")
        return datasets

    return new_datasets

def excludeSamples(excludes, datasets):
    if excludes:
        if isinstance(excludes, list):
            # remove selected datasets
            return list(filter(lambda x: x not in selectSamples(excludes, datasets), datasets))
        elif isinstance(excludes, str):
            # remove selected datasets
            return list(filter(lambda x: x not in selectSample(excludes, datasets), datasets))
        else:
            return list(filter(excludes, datasets))
    else:
        return datasets

