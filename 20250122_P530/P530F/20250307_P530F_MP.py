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

sample_name = 'P530F'
wafername = 'P530'
Lsample = 6 #mm
# angle_of_incidence = 45
sample_thick = 0.4215 #mm
Nbounces = Lsample/sample_thick
# epi_path = np.sqrt(2)*Nperiods*well_period
Lpath = Lsample*np.sqrt(2)

numax = 3100
numin = 1200

blues = ['navy','blue','cornflowerblue']
reds = ['darkviolet','fuchsia','crimson']

#adjust with well thicknesses based on Lodo runsheet

samp_meas = MultipassMeas(samp=sample_name)
bg_meas = MultipassMeas(samp='bg_ap_20_gain_4')

base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/P530F/20250307_P530F_MP'

tm_file = os.path.join(base_dir, 'P530F_TM' + '.CSV')
te_file = os.path.join(base_dir, 'P530F_TE' + '.CSV')

te_bg_file = os.path.join(base_dir,'bg_TE_ap_20_gain_4'+'.CSV')
tm_bg_file = os.path.join(base_dir,'bg_TM_ap_20_gain_4'+'.CSV')

_, tm_wavenum, tm_single_beam, _ = load_data(tm_file)
_, te_wavenum, te_single_beam, _ = load_data(te_file)

_, tm_bg_wavenum, tm_bg_single_beam, _ = load_data(tm_bg_file)
_, te_bg_wavenum, te_bg_single_beam, _ = load_data(te_bg_file)

#try using backgrounds with no sample in the path
samp_meas.TM_wavenum = tm_wavenum
samp_meas.TE_wavenum = te_wavenum
samp_meas.TM_single_beam = tm_single_beam
samp_meas.TE_single_beam = te_single_beam

bg_meas.TM_wavenum = tm_bg_wavenum
bg_meas.TE_wavenum = te_bg_wavenum
bg_meas.TM_single_beam = tm_bg_single_beam
bg_meas.TE_single_beam = te_bg_single_beam

longer_samp_meas = MultipassMeas(samp=sample_name)
longer_samp_name = "P530-A"
l_longer = 10 #mm
Nbounces_longer = l_longer/sample_thick
l_path_longer = l_longer*np.sqrt(2)
numin_longer = 1700
numax_longer = numax

base_dir_longer = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250124_P530-A-MP'

tm_file_longer = os.path.join(base_dir_longer, 'P530-MP-P0deg' + '.CSV')
te_file_longer = os.path.join(base_dir_longer, 'P530-A-MP-P90deg' + '.CSV')

_, tm_wavenum, tm_single_beam, _ = load_data(tm_file_longer)
_, te_wavenum, te_single_beam, _ = load_data(te_file_longer)

#try using backgrounds with no sample in the path
longer_samp_meas.TM_wavenum = tm_wavenum
longer_samp_meas.TE_wavenum = te_wavenum
longer_samp_meas.TM_single_beam = tm_single_beam
longer_samp_meas.TE_single_beam = te_single_beam

# SP_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250122_P530-A-SP'
# SP_name = 'SP 230 deg'
# SP_meas = MultipassMeas(samp=SP_name)
# #
# tm_sp_file = os.path.join(SP_dir, 'P530-SP-rot230deg-P0deg' + '.CSV')
# te_sp_file = os.path.join(SP_dir, 'P530-SP-rot230deg-P90deg' + '.CSV')
#
# _, tm_sp_wavenum, tm_sp_single_beam, _ = load_data(tm_sp_file)
# _, te_sp_wavenum, te_sp_single_beam, _ = load_data(te_sp_file)
#
# SP_meas.TE_single_beam=te_sp_single_beam
# SP_meas.TM_single_beam=tm_sp_single_beam
# SP_meas.TE_wavenum=te_sp_wavenum
# SP_meas.TM_wavenum=tm_sp_wavenum

fig, axs = plt.subplots(1, 3, figsize=(14, 8))

axs[0].plot(samp_meas.TE_wavenum, samp_meas.TE_single_beam, label= sample_name + 'TE MP',color=blues[0])
axs[0].plot(samp_meas.TM_wavenum , samp_meas.TM_single_beam , label= sample_name + 'TM MP', color=reds[0])

axs[0].plot(longer_samp_meas.TE_wavenum, longer_samp_meas.TE_single_beam, label= longer_samp_name + " TE MP",color=blues[1])
axs[0].plot(longer_samp_meas.TM_wavenum , longer_samp_meas.TM_single_beam , label= longer_samp_name+" TM MP", color=reds[1])

axs[0].plot(bg_meas.TE_wavenum, bg_meas.TE_single_beam, label= 'TE bg',color=blues[2])
axs[0].plot(bg_meas.TM_wavenum , bg_meas.TM_single_beam , label= 'TM bg', color=reds[2])

# axs[0].plot(SP_meas.TE_wavenum, SP_meas.TE_single_beam, label= 'TE ' + SP_name,color=blues[1])
# axs[0].plot(SP_meas.TM_wavenum , SP_meas.TM_single_beam , label= 'TM ' + SP_name,color=reds[1])


# axs[1].plot(samp_meas.TE_wavenum, normalize(samp_meas.TE_single_beam), label= 'TE',color=blues[0])
# axs[1].plot(samp_meas.TM_wavenum , normalize(samp_meas.TM_single_beam) , label= 'TM',color=reds[0])

# axs[1].plot(SP_meas.TE_wavenum, normalize(SP_meas.TE_single_beam), label= 'TE ' + SP_name, color=blues[1])
# axs[1].plot(SP_meas.TM_wavenum , normalize(SP_meas.TM_single_beam) , label= 'TM ' + SP_name,  color=reds[1])

# axs[0].plot(bg_meas.TE_wavenum, samp_meas.TE_reshaped, label= 'TE  reshaped',
#                         color='dodgerblue')
# axs[0].plot(bg_meas.TM_wavenum , samp_meas.TM_reshaped , label= 'TM  reshaped',
#                         color='sienna')

#now do the masking

mask_samp = (samp_meas.TE_wavenum > numin) & (samp_meas.TE_wavenum < numax)
#
mask_samp_longer = (longer_samp_meas.TE_wavenum > numin_longer) & (longer_samp_meas.TE_wavenum < numax_longer)
# # Apply the mask to the array
#
mask_bg = (bg_meas.TE_wavenum > numin) & (bg_meas.TE_wavenum < numax)

#average the signal to compare


axs[0].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[0].set_ylabel("Single Beam",fontsize=12)
axs[1].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[1].set_ylabel("TM/TE",fontsize=12)
axs[2].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[2].set_ylabel(r"$\alpha_{ISB} \times L_{path}=-\ln (\frac{I_{out,TM}}{I_{out,TE}}) + \ln(\frac{I_{bg,TM}}{I_{bg,TE}})$",fontsize=12)
# theta_variation_title =sample_name + ' compared to ' + SP_name
theta_variation_title = sample_name + ' Multipass'
axs[0].set_title(theta_variation_title)
axs[0].grid()
transmission_ratio_title = "SNR mask " + str(numin) + r"$ < \nu < $" + str(numax) + r" ${cm}^{-1}$"
axs[1].set_title(transmission_ratio_title)
axs[1].grid()
axs[2].set_title("Absorption")
axs[2].grid()
#
fit_plot_title = theta_variation_title + ' norms'
# fit_plot_title = sample_name+ " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax) + r" ${cm}^{-1}$"

# axs[1].set_title(fit_plot_title)
# axs[1].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
# axs[1].set_ylabel('Normalized by peak absorption',fontsize=12)
# theta_variation_title_polarization_ratios = theta_variation_title + ' polarization ratios'
#
samp_meas.TE_masked = samp_meas.TE_single_beam[mask_samp]
samp_meas.TM_masked = samp_meas.TM_single_beam[mask_samp]
samp_meas.TM_wavenum_masked = samp_meas.TM_wavenum[mask_samp]
samp_meas.TE_wavenum_masked = samp_meas.TE_wavenum[mask_samp]

longer_samp_meas.TE_masked = longer_samp_meas.TE_single_beam[mask_samp_longer]
longer_samp_meas.TM_masked = longer_samp_meas.TM_single_beam[mask_samp_longer]
longer_samp_meas.TM_wavenum_masked = longer_samp_meas.TM_wavenum[mask_samp_longer]
longer_samp_meas.TE_wavenum_masked = longer_samp_meas.TE_wavenum[mask_samp_longer]

bg_meas.TE_masked = bg_meas.TE_single_beam[mask_bg]
bg_meas.TM_masked = bg_meas.TM_single_beam[mask_bg]
bg_meas.TM_wavenum_masked = bg_meas.TM_wavenum[mask_bg]
bg_meas.TE_wavenum_masked = bg_meas.TE_wavenum[mask_bg]
#
path_len_label = r"$N_{bounces} = %0.2f $" % (Nbounces) + r"$,L_{path} = %0.2f mm$" % (Lpath) + r" masked"

axs[1].plot(samp_meas.TE_wavenum_masked, samp_meas.TM_masked/samp_meas.TE_masked, label=path_len_label,
                        color='green')
axs[1].plot(bg_meas.TE_wavenum_masked, bg_meas.TM_masked/bg_meas.TE_masked, label= 'no sample masked ',
                        color='y')
# axs[1].plot(longer_samp_meas.TE_wavenum_masked, longer_samp_meas.TM_masked/longer_samp_meas.TE_masked, label= longer_samp_name + "TM/TE SNR mask " + str(numin_longer) + r"$ < \nu < $" + str(numax_longer) + r" ${cm}^{-1}$",
#                         color='aqua')
path_len_label_longer = r"$N_{bounces} = %0.2f $" % (Nbounces_longer) + r"$,L_{path} = %0.2f mm$" % (l_path_longer) + r" masked"
axs[1].plot(longer_samp_meas.TE_wavenum_masked, longer_samp_meas.TM_masked/longer_samp_meas.TE_masked, label= path_len_label_longer,
                        color='aqua')

offset = np.log(bg_meas.TM_masked / bg_meas.TE_masked)

alpha_ISB = -np.log(samp_meas.TM_masked / samp_meas.TE_masked) + offset

axs[2].plot(samp_meas.TE_wavenum_masked, alpha_ISB, label=sample_name,
              color=reds[0])

axs[0].legend()
axs[0].legend(prop={"size":14})
axs[1].legend()
axs[1].legend(prop={"size":10})
axs[2].legend()
axs[2].legend(prop={"size":14})

plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + longer_samp_name + 'raw data and absorption' + '.svg')
plt.savefig(save_title)
plt.show()