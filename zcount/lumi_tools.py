import csv
import json
import hist
import numpy as np

from utilities import logging

logger = logging.child_logger(__name__)


def get_run_lumi_axes(filename_lumi, filename_json=None):
    # make histogram for run & lumisection with luminosity as bin content; each lumisection comes once, use lumi as weight
    # return the min and max run number, and max lumi section found in the lumi file (that is also present in the json file)
    runs = []
    lumis = []

    with open(filename_lumi) as lumicsv:
        reader = csv.reader(lumicsv)
        for row in reader:
            if row[0][0]=="#":
                continue

            run, _ = row[0].split(":")
            lumi, _ = row[1].split(":")
            
            run = int(run)
            lumi = int(lumi)
            
            runs.append(run)
            lumis.append(lumi)
    logger.info(f"Found runs {min(runs)} to {max(runs)} with lumis up to {max(lumis)}")

    if filename_json is not None:
        logger.info("filter json")
        with open(filename_json) as jsonfile:
            jsondata = json.load(jsonfile)
        runs_json = []
        lumis_json = []
        for run,lumipairs in jsondata.items():
            for lumipair in lumipairs:
                runs_json.append(int(run))
                lumis_json.append(int(lumipair[1]))

        lumis = [lumis_json[runs_json.index(r)] for r,l in zip(runs,lumis) if r in runs_json]
        runs = [r for r in runs if r in runs_json]

        logger.info(f"Found runs {min(runs)} to {max(runs)} with lumis up to {max(lumis)}")

    lumi_max = max(lumis)
    runs = [r for r in set(runs)]
    runs = np.array(runs)
    runs.sort()
    runs = np.append(runs, runs[-1]+1)

    axis_run = hist.axis.Variable(runs, underflow=False, overflow=False, name = "run")
    axis_lumi = hist.axis.Regular(lumi_max, 0, lumi_max+1, underflow=False, overflow=False, name = "lumi")

    return axis_run, axis_lumi