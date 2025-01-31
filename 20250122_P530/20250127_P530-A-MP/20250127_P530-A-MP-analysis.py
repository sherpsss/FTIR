import pandas as pd
import scipy.optimize as opt

import os
# import scipy
import matplotlib.pyplot as plt
import numpy as np
from FTIR_analysis_helpers import MultipassMeas
from FTIR_analysis_helpers import load_data
from FTIR_analysis_helpers import fitLorentzPlot

def normalize(arr):
    amp = np.max(arr)-np.min(arr)
    yzerod = arr-np.min(arr)
    return yzerod/amp

sample_name = 'P530-A-MP'
wafername = 'P530'
Lsample = 10 #mm
# angle_of_incidence = 45
sample_thick = 0.5 #mm
Nbounces = Lsample/sample_thick
well_period = 320e-5 #mm
Nperiods = 29
epi_path = np.sqrt(2)*Nperiods*well_period
Lpath = epi_path*Nbounces

numax = 3100
numin = 750

blues = ['cornflowerblue','mediumblue','blue']
reds = ['crimson','m','magenta']

#adjust with well thicknesses based on Lodo runsheet

samp_meas = MultipassMeas(samp=sample_name)

base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250124_P530-A-MP'

tm_file = os.path.join(base_dir, 'P530-MP-P0deg' + '.CSV')
te_file = os.path.join(base_dir, 'P530-A-MP-P90deg' + '.CSV')

_, tm_wavenum, tm_single_beam, _ = load_data(tm_file)
_, te_wavenum, te_single_beam, _ = load_data(te_file)

#try using backgrounds with no sample in the path
samp_meas.TM_wavenum = tm_wavenum
samp_meas.TE_wavenum = te_wavenum
samp_meas.TM_single_beam = tm_single_beam
samp_meas.TE_single_beam = te_single_beam

SP_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250122_P530-A-SP'
SP_name = 'SP 230 deg'
SP_meas = MultipassMeas(samp=SP_name)
#
tm_sp_file = os.path.join(SP_dir, 'P530-SP-rot230deg-P0deg' + '.CSV')
te_sp_file = os.path.join(SP_dir, 'P530-SP-rot230deg-P90deg' + '.CSV')

_, tm_sp_wavenum, tm_sp_single_beam, _ = load_data(tm_sp_file)
_, te_sp_wavenum, te_sp_single_beam, _ = load_data(te_sp_file)

SP_meas.TE_single_beam=te_sp_single_beam
SP_meas.TM_single_beam=tm_sp_single_beam
SP_meas.TE_wavenum=te_sp_wavenum
SP_meas.TM_wavenum=tm_sp_wavenum

fig, axs = plt.subplots(1, 2, figsize=(14, 8))

axs[0].plot(samp_meas.TE_wavenum, samp_meas.TE_single_beam, label= 'TE MP',color=blues[0])
# axs[0].plot(samp_meas.TM_wavenum , samp_meas.TM_single_beam , label= 'TM MP', color=reds[0])

axs[0].plot(SP_meas.TE_wavenum, SP_meas.TE_single_beam, label= 'TE ' + SP_name,color=blues[1])
# axs[0].plot(SP_meas.TM_wavenum , SP_meas.TM_single_beam , label= 'TM ' + SP_name,color=reds[1])


axs[1].plot(samp_meas.TE_wavenum, normalize(samp_meas.TE_single_beam), label= 'TE',color=blues[0])
# axs[1].plot(samp_meas.TM_wavenum , normalize(samp_meas.TM_single_beam) , label= 'TM',color=reds[0])

axs[1].plot(SP_meas.TE_wavenum, normalize(SP_meas.TE_single_beam), label= 'TE ' + SP_name, color=blues[1])
# axs[1].plot(SP_meas.TM_wavenum , normalize(SP_meas.TM_single_beam) , label= 'TM ' + SP_name,  color=reds[1])

# axs[0].plot(bg_meas.TE_wavenum, samp_meas.TE_reshaped, label= 'TE  reshaped',
#                         color='dodgerblue')
# axs[0].plot(bg_meas.TM_wavenum , samp_meas.TM_reshaped , label= 'TM  reshaped',
#                         color='sienna')

#now do the masking

# mask_samp = (samp_meas.TE_wavenum_reshaped > numin) & (samp_meas.TE_wavenum_reshaped < numax)
# #
# # # Apply the mask to the array
# #
# mask_bg = (bg_meas.TE_wavenum > numin) & (bg_meas.TE_wavenum < numax)



#average the signal to compare


axs[0].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[0].set_ylabel("Single Beam",fontsize=12)
# theta_variation_title =sample_name + ' compared to ' + SP_name
theta_variation_title = wafername + ' multipass vs. single pass'
axs[0].set_title(theta_variation_title)
#
fit_plot_title = theta_variation_title + ' norms'
# fit_plot_title = sample_name+ " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax) + r" ${cm}^{-1}$"

axs[1].set_title(fit_plot_title)
axs[1].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[1].set_ylabel('Normalized by peak absorption',fontsize=12)
# theta_variation_title_polarization_ratios = theta_variation_title + ' polarization ratios'
#
# samp_meas.TE_masked = samp_meas.TE_reshaped[mask_samp]
# samp_meas.TM_masked = samp_meas.TM_reshaped[mask_samp]
# samp_meas.TM_wavenum_masked = samp_meas.TM_wavenum_reshaped[mask_samp]
# samp_meas.TE_wavenum_masked = samp_meas.TE_wavenum_reshaped[mask_samp]
#
# axs[1].plot(samp_meas.TE_wavenum_masked, samp_meas.TM_masked/samp_meas.TE_masked, label= 'TM/TE with sample masked ',
#                         color='green')
# axs[1].plot(bg_meas.TE_wavenum[mask_bg], bg_meas.TM_single_beam[mask_bg]/bg_meas.TE_single_beam[mask_bg], label= 'TM/TE no sample masked ',
#                         color='y')

axs[0].legend()
axs[0].legend(prop={"size":14})

plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + ' TE raw data' + '.svg')
plt.savefig(save_title)
plt.show()
