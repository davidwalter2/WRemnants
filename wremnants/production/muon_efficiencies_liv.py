import ROOT

import narf
from wremnants.utilities import common
from wums import logging

logger = logging.child_logger(__name__)

narf.clingutils.Declare('#include "muon_reco_efficiency.hpp"')

data_dir = common.data_dir


def make_muon_reco_efficiency_helper(
    era=None,
):

    filename = data_dir + "muonSF/allSmooth_GtoHout_vtxAgnIso_altBkg.root"
    dataname = f"SF_nomiAndAlt_GtoH_reco_plus_altBkg"

    fdata = ROOT.TFile.Open(filename)
    datahist = fdata.Get(dataname)

    datahist.SetDirectory(0)  ### i don't know what this does
    fdata.Close()

    logger.debug(f"Reading SF file-{filename}")
    datahist.SetName(f"reco_SF_{era}")
    datahist.SetTitle("")
    logger.debug("")
    logger.debug(f"reco scale factors for era {era}")
    logger.debug(
        [
            datahist.GetBinContent(i, j)
            for i in range(1, datahist.GetNbinsX() + 1)
            for j in range(1, datahist.GetNbinsY() + 1)
        ]
    )
    logger.debug("")

    helper = ROOT.wrem.muon_reco_efficiency_helper(datahist)

    return helper


def make_muon_tracking_efficiency_helper(
    era=None,
):

    filename = data_dir + "muonSF/allSmooth_GtoHout_vtxAgnIso_altBkg.root"
    dataname = f"SF_nomiAndAlt_GtoH_tracking_plus_altBkg"

    fdata = ROOT.TFile.Open(filename)
    datahist = fdata.Get(dataname)

    datahist.SetDirectory(0)  ### i don't know what this does
    fdata.Close()

    logger.debug(f"Reading SF file-{filename}")

    datahist.SetName(
        f"reco_SF_{era}"
    )  ### i do want to keep this because tehcnically i should be applying this by ero
    datahist.SetTitle("")
    logger.debug("")
    logger.debug(f"reco scale factors for era {era}")
    logger.debug(
        [
            datahist.GetBinContent(i, j)
            for i in range(1, datahist.GetNbinsX() + 1)
            for j in range(1, datahist.GetNbinsY() + 1)
        ]
    )
    logger.debug("")

    ### the problem is that i need to be passing two values in here.
    ### but in other ones they seem to pass a histogram in here, not specific values.
    helper = ROOT.wrem.muon_reco_efficiency_helper(datahist)

    return helper
