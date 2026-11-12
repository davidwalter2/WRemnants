import matplotlib.pyplot as plt
import h5py
from work.WRemnants.utilities.io_tools import input_tools
from work.WRemnants.wums.wums.boostHistHelpers import broadcastSystHist, multiplyHists, divideHists
import numpy as np
from scipy.optimize import curve_fit


file_in = '/work/submit/jbenke/WRemnants/scripts/histmakers/'
# file_out = '/home/submit/jbenke/public_html/LIV_mass_evo/'
file_out = '/home/submit/jbenke/public_html/'


def make_plot(data_all, plotname, legend_all = ["MC", "Data"], ylim = []):
    plt.clf()
    for i in range(len(legend_all)):
        data_all[i].plot1d()
    plt.xlim([0.01, 24])
    if len(ylim) != 0:
        plt.ylim(ylim)
    # if xaxis != '':
    #     plt.xlabel(xaxis)
    plt.title(plotname)
    plt.legend(legend_all)
    
    plt.savefig(file_out + plotname + '.png')
file_in_name = file_in + 'mz_dilepton_liv_scetlib_dyturboCorr_maxFiles_20.hdf5'
h5file = h5py.File(file_in_name, 'r')
results = input_tools.load_results_h5py(h5file)





time_mll_data = results['dataPostVFP']['output']['time_mll'].get()
mll_data = time_mll_data.project('mll')
time_data = time_mll_data.project('time')

mass_mc = results['ZmumuPostVFP']['output']['mll'].get() 
# mass_mc = results['ZmumuPostVFP']['output']['pass_reco_pass_gen'].get().project("mll") + results['ZmumuPostVFP']['output']['pass_reco_fail_gen'].get().project("mll")

lumi_scaling = results['dataPostVFP']['lumi_outout']['time'].get()

weightsum = results['ZmumuPostVFP']['weight_sum']
cross_sec = results["ZmumuPostVFP"]["dataset"]["xsec"]

def mc_corrections(mc_results):
    mc_results /= weightsum
    mc_results *= cross_sec
    mc_results *= 1000
    return mc_results

stight_strig_data = results['dataPostVFP']['output']['time_mll_stight_strig'].get()
dtight_strig_data = results['dataPostVFP']['output']['time_mll_dtight_strig'].get()

stight_strig_mc = results['ZmumuPostVFP']['output']['mll_stight_strig'].get()
dtight_strig_mc = results['ZmumuPostVFP']['output']['mll_dtight_strig'].get()

generator = results['ZmumuPostVFP']['output']['pass_gen'].get()

mass_mc = mc_corrections(mass_mc)
stight_strig_mc = mc_corrections(stight_strig_mc)
dtight_strig_mc = mc_corrections(dtight_strig_mc)
generator = mc_corrections(generator)

mass_mc_2d = broadcastSystHist(mass_mc, time_mll_data)
stight_strig_mc_2d = broadcastSystHist(stight_strig_mc, time_mll_data)
dtight_strig_mc_2d = broadcastSystHist(dtight_strig_mc, time_mll_data)
generator_2d = broadcastSystHist(generator, time_mll_data)

mass_mc_2d = multiplyHists(mass_mc_2d, lumi_scaling)
stight_strig_mc_2d = multiplyHists(stight_strig_mc_2d, lumi_scaling)
dtight_strig_mc_2d = multiplyHists(dtight_strig_mc_2d, lumi_scaling)
generator_2d =  multiplyHists(generator_2d, lumi_scaling)


for i in range(len(time_data.values())):
    ### UGH THIS DOESN'T WORK ANYMORE
    time_mc_proj = mass_mc_2d[:, i]
    time_data_proj = time_mll_data[:, i]
    generator_proj = generator_2d[:, i]
    
    stight_strig_mc_proj = stight_strig_mc_2d[:, i]
    dtight_strig_mc_proj = dtight_strig_mc_2d[:, i]
    
    stight_strig_data_proj = stight_strig_data[:, i]
    dtight_strig_data_proj = dtight_strig_data[:, i]
    
    if i == 0: 
        mass_sum_mc = time_mc_proj
        mass_sum_data = time_data_proj
        mass_sum_stight_strig_mc = stight_strig_mc_proj
        mass_sum_dtight_strig_mc = dtight_strig_mc_proj
        mass_sum_generator = generator_proj

        mass_sum_stight_strig_data = stight_strig_data_proj
        mass_sum_dtight_strig_data = dtight_strig_data_proj
    else:
        mass_sum_mc += time_mc_proj
        mass_sum_data += time_data_proj
        mass_sum_stight_strig_mc += stight_strig_mc_proj
        mass_sum_dtight_strig_mc += dtight_strig_mc_proj

        mass_sum_stight_strig_data += stight_strig_data_proj
        mass_sum_dtight_strig_data += dtight_strig_data_proj
        mass_sum_generator += generator_proj

    
    time_mc_proj.plot1d(label = 'MC: 2 tightID, 2 trig')
    time_data_proj.plot1d(label = 'Data: 2 tightID, 2 trig')
    # stight_strig_mc_proj.plot1d(label = ' MC: 1 tightID, 1 trig')
    # stight_strig_data_proj.plot1d(label = 'Data: 1 tightID, 1 trig')
    # dtight_strig_mc_proj.plot1d(label = 'MC: 2 tightID, 1 trig')
    # dtight_strig_data_proj.plot1d(label = 'Data: 2 tightID, 1 trig')
    # generator_proj.plot1d("label = Generator")
    
    make_plot([time_mc_proj, time_data_proj], "mass_bin_" + str(i))

make_plot([mass_sum_mc, mass_sum_data, mass_sum_stight_strig_mc,  mass_sum_stight_strig_data,mass_sum_dtight_strig_mc, mass_sum_dtight_strig_data, mass_sum_generator], "mass_sum_all", legend_all = ['MC: 2 tightID, 2 trig', 'Data: 2 tightID, 2 trig', 'MC: 1 tightID, 1 trig', 'Data: 1 tightID, 1 trig', 'MC: 2 tightID, 1 trig', 'Data: 2 tightID, 1 trig', 'Generator'])

make_plot([mass_sum_generator], "generator", legend_all = ["Generator"])

make_plot([mass_sum_mc, mass_sum_data], "mass_sum_new")
ratio_mass_sum = divideHists(mass_sum_mc, mass_sum_data)
make_plot([ratio_mass_sum], "mass_sum_ratio", ["MC/Data"], [1, 1.5])














###### this worked the morning of 7/2 ####    
# data_file =
'''
time_mll_data = results['dataPostVFP']['output']['time_mll'].get()
mll_data = time_mll_data.project('mll')
time_data = time_mll_data.project('time')

mass_mc = results['ZmumuPostVFP']['output']['mll'].get() 
lumi_scaling = results['dataPostVFP']['lumi_outout']['time'].get()

weightsum = results['ZmumuPostVFP']['weight_sum']
cross_sec = results["ZmumuPostVFP"]["dataset"]["xsec"]

def mc_corrections(mc_results):
    mc_results /= weightsum
    mc_results *= cross_sec
    mc_results *= 1000
    return mc_results

stight_strig_data = results['dataPostVFP']['output']['time_mll_stight_strig'].get() 
dtight_strig_data = results['dataPostVFP']['output']['time_mll_dtight_strig'].get() 

stight_strig_mc = results['ZmumuPostVFP']['output']['mll_stight_strig'].get() 
dtight_strig_mc = results['ZmumuPostVFP']['output']['mll_dtight_strig'].get() 

generator = results['ZmumuPostVFP']['output']['gen_mll'].get()

mass_mc = mc_corrections(mass_mc)
stight_strig_mc = mc_corrections(stight_strig_mc)
dtight_strig_mc = mc_corrections(dtight_strig_mc)
generator = mc_corrections(generator)

mass_mc_2d = broadcastSystHist(mass_mc, time_mll_data)
stight_strig_mc_2d = broadcastSystHist(stight_strig_mc, time_mll_data)
dtight_strig_mc_2d = broadcastSystHist(dtight_strig_mc, time_mll_data)
generator_2d = broadcastSystHist(generator, time_mll_data)

mass_mc_2d = multiplyHists(mass_mc_2d, lumi_scaling)
stight_strig_mc_2d = multiplyHists(stight_strig_mc_2d, lumi_scaling)
dtight_strig_mc_2d = multiplyHists(dtight_strig_mc_2d, lumi_scaling)
generator_2d =  multiplyHists(generator_2d, lumi_scaling)


for i in range(len(time_data.values())):
    ### UGH THIS DOESN'T WORK ANYMORE
    time_mc_proj = mass_mc_2d[:, i]
    time_data_proj = time_mll_data[:, i]
    generator_proj = generator_2d[:, i]
    
    stight_strig_mc_proj = stight_strig_mc_2d[:, i]
    dtight_strig_mc_proj = dtight_strig_mc_2d[:, i]
    
    stight_strig_data_proj = stight_strig_data[:, i]
    dtight_strig_data_proj = dtight_strig_data[:, i]
    
    if i == 0: 
        mass_sum_mc = time_mc_proj
        mass_sum_data = time_data_proj
        mass_sum_stight_strig_mc = stight_strig_mc_proj
        mass_sum_dtight_strig_mc = dtight_strig_mc_proj
        mass_sum_generator = generator_proj

        mass_sum_stight_strig_data = stight_strig_data_proj
        mass_sum_dtight_strig_data = dtight_strig_data_proj
    else:
        mass_sum_mc += time_mc_proj
        mass_sum_data += time_data_proj
        mass_sum_stight_strig_mc += stight_strig_mc_proj
        mass_sum_dtight_strig_mc += dtight_strig_mc_proj

        mass_sum_stight_strig_data += stight_strig_data_proj
        mass_sum_dtight_strig_data += dtight_strig_data_proj
        mass_sum_generator += generator_proj

    
    time_mc_proj.plot1d(label = 'MC: 2 tightID, 2 trig')
    time_data_proj.plot1d(label = 'Data: 2 tightID, 2 trig')
    # stight_strig_mc_proj.plot1d(label = ' MC: 1 tightID, 1 trig')
    # stight_strig_data_proj.plot1d(label = 'Data: 1 tightID, 1 trig')
    # dtight_strig_mc_proj.plot1d(label = 'MC: 2 tightID, 1 trig')
    # dtight_strig_data_proj.plot1d(label = 'Data: 2 tightID, 1 trig')
    # generator_proj.plot1d("label = Generator")
    
    make_plot([time_mc_proj, time_data_proj], "mass_bin_" + str(i))

make_plot([mass_sum_mc, mass_sum_data, mass_sum_stight_strig_mc,  mass_sum_stight_strig_data,mass_sum_dtight_strig_mc, mass_sum_dtight_strig_data, mass_sum_generator], "mass_sum_all", legend_all = ['MC: 2 tightID, 2 trig', 'Data: 2 tightID, 2 trig', 'MC: 1 tightID, 1 trig', 'Data: 1 tightID, 1 trig', 'MC: 2 tightID, 1 trig', 'Data: 2 tightID, 1 trig', 'Generator'])

make_plot([mass_sum_generator], "generator", legend_all = ["Generator"])

make_plot([mass_sum_mc, mass_sum_data], "mass_sum_new")
ratio_mass_sum = divideHists(mass_sum_mc, mass_sum_data)
make_plot([ratio_mass_sum], "mass_sum_ratio", ["MC/Data"], [1, 1.5])


'''