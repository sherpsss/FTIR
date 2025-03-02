import pandas as pd
import scipy.optimize as opt

import os
# import scipy
import matplotlib.pyplot as plt
import numpy as np
from FTIR_analysis_helpers import MultipassMeas
from FTIR_analysis_helpers import load_data
from FTIR_analysis_helpers import fitLorentzPlot

sample_name = 'P12-3-24-1'
Lsample = 10 #mm
# angle_of_incidence = 45
sample_thick = 0.5 #mm
Nbounces = Lsample/sample_thick
well_period = 320e-5 #mm
Nperiods = 29
epi_path = np.sqrt(2)*Nperiods*well_period
Lpath = epi_path*Nbounces

numax = 1713
numin = 750


#adjust with well thicknesses based on Lodo runsheet

samp_meas = MultipassMeas(samp=sample_name)

base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/P12-3-24-1'

tm_file = os.path.join(base_dir, 'P0deg' + '.CSV')
te_file = os.path.join(base_dir, 'P90deg' + '.CSV')

_, tm_wavenum, tm_single_beam, _ = load_data(tm_file)
_, te_wavenum, te_single_beam, _ = load_data(te_file)

#try using backgrounds with no sample in the path
samp_meas.TM_wavenum = tm_wavenum
samp_meas.TE_wavenum = te_wavenum
samp_meas.TM_single_beam = tm_single_beam
samp_meas.TE_single_beam = te_single_beam

background_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20241211_theta_i_sweep_singlepass_342A/backgrounds'
bg_name = 'no_samp'
bg_meas = MultipassMeas(samp=bg_name)

tm_bg_file = os.path.join(background_dir, 'background_P0deg' + '.CSV')
te_bg_file = os.path.join(background_dir, 'background_P90deg' + '.CSV')

_, tm_bg_wavenum, tm_bg_single_beam, _ = load_data(tm_bg_file)
_, te_bg_wavenum, te_bg_single_beam, _ = load_data(te_bg_file)

bg_meas.TE_single_beam=te_bg_single_beam
bg_meas.TM_single_beam=tm_bg_single_beam
bg_meas.TE_wavenum=te_bg_wavenum
bg_meas.TM_wavenum=tm_bg_wavenum

fig, axs = plt.subplots(1, 2, figsize=(14, 8))

axs[0].plot(samp_meas.TE_wavenum, samp_meas.TE_single_beam, label= 'TE',
                        color='blue')
axs[0].plot(samp_meas.TM_wavenum , samp_meas.TM_single_beam , label= 'TM',
                        color='red')

axs[0].plot(bg_meas.TE_wavenum, bg_meas.TE_single_beam, label= 'TE bg ' + bg_name,
                        color='c')
axs[0].plot(bg_meas.TM_wavenum , bg_meas.TM_single_beam , label= 'TM bg' + bg_name,
                        color='m')

samp_meas.TE_reshaped = samp_meas.TE_single_beam.reshape(-1,2).mean(axis=1)
samp_meas.TM_reshaped = samp_meas.TM_single_beam.reshape(-1,2).mean(axis=1)
samp_meas.TE_wavenum_reshaped = samp_meas.TE_wavenum.reshape(-1,2).mean(axis=1)
samp_meas.TM_wavenum_reshaped = samp_meas.TM_wavenum.reshape(-1,2).mean(axis=1)

# axs[0].plot(bg_meas.TE_wavenum, samp_meas.TE_reshaped, label= 'TE  reshaped',
#                         color='dodgerblue')
# axs[0].plot(bg_meas.TM_wavenum , samp_meas.TM_reshaped , label= 'TM  reshaped',
#                         color='sienna')

#now do the masking

mask_samp = (samp_meas.TE_wavenum_reshaped > numin) & (samp_meas.TE_wavenum_reshaped < numax)
#
# # Apply the mask to the array
#
mask_bg = (bg_meas.TE_wavenum > numin) & (bg_meas.TE_wavenum < numax)



#average the signal to compare

axislabelsfont=30
legendfont=30
ticksize=30
titlesize=34


axs[0].set_xlabel('Wavenumber (cm^-1)',fontsize=axislabelsfont)
axs[0].set_ylabel("Single Beam",fontsize=axislabelsfont)
theta_variation_title =sample_name
axs[0].set_title(theta_variation_title)

fit_plot_title = sample_name+ " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax) + r" ${cm}^{-1}$"
axs[1].set_title(fit_plot_title)
axs[1].set_xlabel(r"$Wavenumber [{cm}^{-1}]$",fontsize=axislabelsfont)
axs[1].set_ylabel('Transmission Ratio',fontsize=axislabelsfont)
theta_variation_title_polarization_ratios = theta_variation_title + ' polarization ratios'

samp_meas.TE_masked = samp_meas.TE_reshaped[mask_samp]
samp_meas.TM_masked = samp_meas.TM_reshaped[mask_samp]
samp_meas.TM_wavenum_masked = samp_meas.TM_wavenum_reshaped[mask_samp]
samp_meas.TE_wavenum_masked = samp_meas.TE_wavenum_reshaped[mask_samp]

axs[1].plot(samp_meas.TE_wavenum_masked, samp_meas.TM_masked/samp_meas.TE_masked, label= 'TM/TE with sample masked ',
                        color='green')
axs[1].plot(bg_meas.TE_wavenum[mask_bg], bg_meas.TM_single_beam[mask_bg]/bg_meas.TE_single_beam[mask_bg], label= 'TM/TE no sample masked ',
                        color='y')

axs[0].legend()
axs[0].legend(prop={"size":14})
axs[1].legend()
axs[1].legend(prop={"size":14})
# plt.legend(["age", "number"], prop = { "size": 20 }, loc ="upper left")
plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'raw data and ratios' + '.svg')
plt.savefig(save_title)
# plt.show()

#fit figure

fig_fits, axs_fits = plt.subplots(figsize=(10, 8))
axs_fits.set_xlabel(r"Wavenumber $[{cm}^{-1}]$",fontsize=axislabelsfont)

offset = np.log(bg_meas.TM_single_beam[mask_bg]/bg_meas.TE_single_beam[mask_bg])

alpha_ISB= (-np.log(samp_meas.TM_masked/samp_meas.TE_masked)+offset)/Nperiods

# axs_fits.plot(samp_meas.TE_wavenum_masked, alpha_ISB, label= r"$\alpha_{ISB}$ per well",
                        # color='green')
axs_fits.plot(samp_meas.TE_wavenum_masked, alpha_ISB, color='green')

# #do fits
# numins = [840,1012,1700,1881]
# numaxs = [1039,1700,1780,1910]
# numins = [860,1012,1700,1881]
# numaxs = [1039,1700,1780,1910]
# kappa_nu_guesses = [35,50,10,5]
numins = [860,1012]
numaxs = [1039,1700]
kappa_nu_guesses = [35,50]
# numins = [840,1012]
# numaxs = [1039,1700]
for i in range(0,len(numins)):
    nu_range = [numins[i],numaxs[i]]
    kappa_nu_guess = kappa_nu_guesses[i]
    fitLorentzParams = fitLorentzPlot(nu_range, kappa_nu_guess,samp_meas.TE_wavenum_masked, alpha_ISB, axs_fits)
    # fitNormalPlot(nu_range,fitLorentzParams[0],fitLorentzParams[1],samp_meas.TE_wavenum_masked,alpha_ISB,axs_fits)


axs_fits.set_ylabel(r"$\alpha_{ISB}$ per well",fontsize=axislabelsfont)
axs_fits.tick_params(axis='x',labelsize=ticksize)
axs_fits.tick_params(axis='y',labelsize=ticksize)
# fit_plot_title = sample_name+ " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax)
fit_plot_title = r"Example Sample Absorption"
axs_fits.set_title(fit_plot_title,fontsize=titlesize)
axs_fits.legend()
axs_fits.legend(prop={"size":legendfont})
plt.figure(fig_fits)
# plt.yticks(fontsize=20)
# plt.xticks(fontisze=20)
plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'Lorentzian fits' + '.svg')
plt.savefig(save_title)
# plt.xlim(numins[2],numaxs[2])
# save_title = os.path.join(base_dir, sample_name + 'Lorentzian fits 1740cm zoomed' + '.svg')
# plt.savefig(save_title)
# plt.xlim(numins[3],numaxs[3])
# save_title = os.path.join(base_dir, sample_name + 'Lorentzian fits 1895cm zoomed' + '.svg')
# plt.savefig(save_title)
plt.figure(fig_fits)
plt.show()
plt.figure(fig)
plt.show()