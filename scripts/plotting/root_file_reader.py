import pdb

import numpy as np
import ROOT

from wremnants.utilities import common

data_dir = common.data_dir
# root -l allSmooth_GtoHout_vtxAgnIso_altBkg.root
era = "2016_PostVFP"
filepath = data_dir + "muonSF/tagAndProbe/2016/"
filename = filepath + "allEfficiencies_2D_idip_minus.root"
# dataname = "effData_nomiAndAlt_GtoH_trigger_plus"
dataname = "SF_nomiAndAlt_GtoH_idip_plus"
# dataname = "effMC_nomiAndAlt_GtoH_idip_plus"

## scale factor: SF2D_nominal
## data: EffData2D
## mc: EffMC2   there are these things for alternate data and im not sure what those mean but im not going to choose that


fdata = ROOT.TFile.Open(filename)
datahist = fdata.Get(dataname).Clone()

datahist.SetDirectory(0)  ### i don't know what this does
fdata.Close()

nx = datahist.GetNbinsX()
ny = datahist.GetNbinsY()
nz = datahist.GetNbinsZ()

# Create a NumPy array with the same shape
data_array = np.zeros((nx, ny, nz))
### axes: eta-pt-ut
for i in range(1, nx + 1):  # eta
    for j in range(1, ny + 1):  # pt
        for k in range(1, nz + 1):  # value?
            data_array[i - 1, j - 1, k - 1] = datahist.GetBinContent(i, j, k)

    ### i need to get the axis spacing as well
pdb.set_trace()
print(data_array)
