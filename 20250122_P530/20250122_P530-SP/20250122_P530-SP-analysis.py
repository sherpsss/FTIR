import os
# import scipy
import matplotlib.pyplot as plt
import numpy as np
from FTIR_analysis_helpers import MultipassMeas
from FTIR_analysis_helpers import load_data
from FTIR_analysis_helpers import fitLorentzPlot
from FTIR_analysis_helpers import SinglePassMeas
# from ...FTIR_analysis_helpers import MultipassMeas
# from ...FTIR_analysis_helpers import fitLorentzPlot

sample_name = 'P530'
Lsample = 10 #mm
# angle_of_incidence = 45
sample_thick = 0.5 #mm
Nbounces = Lsample/sample_thick
well_period = 320e-5 #mm
Nperiods = 29
epi_path = np.sqrt(2)*Nperiods*well_period
Lpath = epi_path*Nbounces

numax = 3100
numin = 785

#adjust with well thicknesses based on Lodo runsheet

base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250122_P530-A-SP'

angles = [230,240,250]
angle0 = 280

#raw data plots
fig, axs = plt.subplots(1, 3, figsize=(14, 8))

axs[0].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[0].set_ylabel("Single Beam",fontsize=12)
theta_variation_title =sample_name + ' Raw data'
axs[0].set_title(theta_variation_title)

fit_plot_title = sample_name+ " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax) + r" ${cm}^{-1}$"
axs[1].set_title(fit_plot_title)
axs[1].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[1].set_ylabel('Transmission Ratio',fontsize=12)

fit_plot_title = sample_name+ " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax) + r" ${cm}^{-1} Samp vs. bg$"
axs[2].set_title(fit_plot_title)
axs[2].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[2].set_ylabel('Raw Transmission Ratios',fontsize=12)

tm_bg_file = os.path.join(base_dir, 'no_samp_P0deg' + '.CSV')
te_bg_file = os.path.join(base_dir, 'no_samp_P90deg' + '.CSV')

_, tm_bg_wavenum, tm_bg_single_beam, _ = load_data(tm_bg_file)
_, te_bg_wavenum, te_bg_single_beam, _ = load_data(te_bg_file)

bg_meas = SinglePassMeas(ident='background')

bg_meas.TE_single_beam=te_bg_single_beam
bg_meas.TM_single_beam=tm_bg_single_beam
bg_meas.TE_wavenum=te_bg_wavenum
bg_meas.TM_wavenum=tm_bg_wavenum

axs[0].plot(bg_meas.TE_wavenum, bg_meas.TE_single_beam, label='TE bg ',
            color='yellowgreen')
axs[0].plot(bg_meas.TM_wavenum, bg_meas.TM_single_beam, label='TM bg',
            color='darkgoldenrod')

mask_samp = (bg_meas.TE_wavenum > numin) & (bg_meas.TE_wavenum < numax)

bg_meas.TE_masked = bg_meas.TE_single_beam[mask_samp]
bg_meas.TM_masked = bg_meas.TM_single_beam[mask_samp]
bg_meas.TM_wavenum_masked = bg_meas.TM_wavenum[mask_samp]
bg_meas.TE_wavenum_masked = bg_meas.TE_wavenum[mask_samp]

axs[1].plot(bg_meas.TM_wavenum_masked, bg_meas.TM_masked/bg_meas.TE_masked,
            label='TM/TE no sample', color='y')

axislabelsfont=30
legendfont=30
ticksize=30
titlesize=34

#absorption fit plot
fig_fits, axs_fits = plt.subplots(figsize=(12, 8))
axs_fits.set_xlabel(r"Wavenumber $[{cm}^{-1}]$",fontsize=axislabelsfont)
# axs_fits.set_ylabel(r"$\alpha_{ISB} \times L_{path}=-\ln (\frac{I_{out,TM}}{I_{out,TE}}) + \ln(\frac{I_{bg,TM}}{I_{bg,TE}})$",fontsize=12)
axs_fits.set_ylabel(r"$\alpha_{ISB} \times L_{path}$",fontsize=axislabelsfont)
# fit_plot_title = sample_name + " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax)
fit_plot_title = "Complex Multi-Quantum Well Sample Absorption"

axs_fits.set_title(fit_plot_title,fontsize=titlesize)
axs_fits.tick_params(axis='x',labelsize=ticksize)
axs_fits.tick_params(axis='y',labelsize=ticksize)

blues = ['darkblue','mediumblue','blue','cornflowerblue']
reds = ['mediumvioletred','deeppink','hotpink','pink']

numins = [1042,1042,1042]
numaxs = [1272,1272,1272]
nuplot_range =[1000,1272]
kappa_nu_guesses = [30]

for i in range(0,len(angles)):
    angle = angles[i]
    angle_filename = sample_name + '-SP-rot'+ str(angle)+'deg-'
    tm_file = os.path.join(base_dir, angle_filename+'P0deg' + '.CSV')
    te_file = os.path.join(base_dir, angle_filename+'P90deg' + '.CSV')
    _, tm_wavenum, tm_single_beam, _ = load_data(tm_file)
    _, te_wavenum, te_single_beam, _ = load_data(te_file)
    angle = angle0-angle
    angle_meas = SinglePassMeas(ident=str(angle))
    angle_meas.TM_wavenum = tm_wavenum
    angle_meas.TE_wavenum = te_wavenum
    angle_meas.TM_single_beam = tm_single_beam
    angle_meas.TE_single_beam = te_single_beam

    #add the raw data plot
    axs[0].plot(angle_meas.TE_wavenum, angle_meas.TE_single_beam, label='TE '+str(angle) + '$\degree$',
                color=blues[i])
    axs[0].plot(angle_meas.TM_wavenum, angle_meas.TM_single_beam, label='TM '+str(angle) + '$\degree$',
                color=reds[i])

    #do the masking
    angle_meas.TE_masked = angle_meas.TE_single_beam[mask_samp]
    angle_meas.TM_masked = angle_meas.TM_single_beam[mask_samp]
    angle_meas.TM_wavenum_masked = angle_meas.TM_wavenum[mask_samp]
    angle_meas.TE_wavenum_masked = angle_meas.TE_wavenum[mask_samp]

    axs[1].plot(angle_meas.TM_wavenum_masked, angle_meas.TM_masked/angle_meas.TE_masked, label= 'TM/TE '+str(angle) + '$\degree$',
                            color=reds[i])

    axs[2].plot(angle_meas.TM_wavenum_masked, angle_meas.TM_masked/bg_meas.TM_masked, label= 'samp TM '+str(angle) + '$\degree $ / TM background' ,
                            color=reds[i])
    axs[2].plot(angle_meas.TE_wavenum_masked, angle_meas.TE_masked/bg_meas.TE_masked, label= 'samp TE '+str(angle) + '$\degree $ / TE background' ,
                            color=blues[i])

    #calculate the absorption coefficient times path length

    offset = np.log(bg_meas.TM_masked/bg_meas.TE_masked)

    alpha_ISB= -np.log(angle_meas.TM_masked/angle_meas.TE_masked)+offset
    axs_fits.plot(angle_meas.TE_wavenum_masked, alpha_ISB, label= str(angle) + '$\degree$',
                            color=reds[i])

    # do the fit
    nu_range = [numins[i],numaxs[i]]
    kappa_nu_guess = kappa_nu_guesses[0]
    # fitLorentzParams = fitLorentzPlot(nu_range, kappa_nu_guess,angle_meas.TE_wavenum_masked, alpha_ISB, axs_fits,nu_fit_plot_range=nuplot_range)

#get the 0 deg one in there
special_filename = 'P530-SP-top-chip'
special_file = os.path.join(base_dir, special_filename + '.CSV')
_, te_wavenum, te_single_beam, _ = load_data(special_file)
angle = 0
angle_meas = SinglePassMeas(ident=str(angle))
angle_meas.TE_wavenum = te_wavenum
angle_meas.TE_single_beam = te_single_beam
axs[0].plot(angle_meas.TE_wavenum, angle_meas.TE_single_beam, label='TE ' + str(angle) + '$\degree$',
            color=blues[-1])

#try using backgrounds with no sample in the path
#
# axs[0].plot(samp_meas.TE_wavenum, samp_meas.TE_single_beam, label= 'TE',
#                         color='blue')
# axs[0].plot(samp_meas.TM_wavenum , samp_meas.TM_single_beam , label= 'TM',
#                         color='red')
#
# axs[0].plot(bg_meas.TE_wavenum, bg_meas.TE_single_beam, label= 'TE bg ' + bg_name,
#                         color='c')
# axs[0].plot(bg_meas.TM_wavenum , bg_meas.TM_single_beam , label= 'TM bg' + bg_name,
#                         color='m')
#
# #now do the masking
#
# mask_samp = (samp_meas.TE_wavenum_reshaped > numin) & (samp_meas.TE_wavenum_reshaped < numax)
# #
# # # Apply the mask to the array
# #
# mask_bg = (bg_meas.TE_wavenum > numin) & (bg_meas.TE_wavenum < numax)
#
# #average the signal to compare
#
# fit_plot_title = sample_name+ " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax) + r" ${cm}^{-1}$"
# axs[1].set_title(fit_plot_title)
# axs[1].set_xlabel('Wavenumber (cm^-1)')
# axs[1].set_ylabel('Transmission Ratio')
# theta_variation_title_polarization_ratios = theta_variation_title + ' polarization ratios'
#
# samp_meas.TE_masked = samp_meas.TE_reshaped[mask_samp]
# samp_meas.TM_masked = samp_meas.TM_reshaped[mask_samp]
# samp_meas.TM_wavenum_masked = samp_meas.TM_wavenum_reshaped[mask_samp]
# samp_meas.TE_wavenum_masked = samp_meas.TE_wavenum_reshaped[mask_samp]
#
# axs[1].plot(samp_meas.TE_wavenum_masked, samp_meas.TM_masked/samp_meas.TE_masked, label= 'TM/TE with sample masked reshaped',
#                         color='green')
# axs[1].plot(bg_meas.TE_wavenum[mask_bg], bg_meas.TM_single_beam[mask_bg]/bg_meas.TE_single_beam[mask_bg], label= 'TM/TE no sample masked reshaped',
#                         color='y')
#
plt.figure(fig)
axs[0].legend()
axs[0].legend(prop={"size":14})
axs[1].legend()
axs[1].legend(prop={"size":14})
axs[2].legend()
axs[2].legend(prop={"size":14})
plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'raw scans and ratios' + '.svg')
plt.savefig(save_title)
#
# #fit figure
#
# offset = np.log(bg_meas.TM_single_beam[mask_bg]/bg_meas.TE_single_beam[mask_bg])
#
# alpha_ISB= -np.log(samp_meas.TM_masked/samp_meas.TE_masked)+offset
#
# axs_fits.plot(samp_meas.TE_wavenum_masked, alpha_ISB, label= r"$-\ln (\frac{I_{out,TM}}{I_{out,TE}}) + \ln(\frac{I_{bg,TM}}{I_{bg,TE}})$",
#                         color='green')
#
# #do fits
# numins = [840,1012,1700,1881]
# numaxs = [1039,1700,1780,1910]
# numins = [860,1012,1700,1881]
# numaxs = [1039,1700,1780,1910]
# kappa_nu_guesses = [35,50,10,5]
# # numins = [840,1012]
# # numaxs = [1039,1700]
# for i in range(0,len(numins)):
#     nu_range = [numins[i],numaxs[i]]
#     kappa_nu_guess = kappa_nu_guesses[i]
#     fitLorentzParams = fitLorentzPlot(nu_range, kappa_nu_guess,samp_meas.TE_wavenum_masked, alpha_ISB, axs_fits)
#     # FTIR_analysis_helpers.fitNormalPlot(nu_range,fitLorentzParams[0],fitLorentzParams[1],samp_meas.TE_wavenum_masked,alpha_ISB,axs_fits)
#
plt.figure(fig_fits)
axs_fits.legend()
axs_fits.legend(prop={"size":legendfont})
plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'alpha_ISBs' + '.svg')
plt.savefig(save_title)

plt.figure(fig)
plt.show()
plt.figure(fig_fits)
plt.show()