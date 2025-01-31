import os
import matplotlib.pyplot as plt
import numpy as np
from FTIR_analysis_helpers import load_data
from FTIR_analysis_helpers import SinglePassMeas


sample_name = '342A'

numax = 3100
numin = 850

#adjust with well thicknesses based on Lodo runsheet

base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20241211_theta_i_sweep_singlepass_342A'
bg_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20241211_theta_i_sweep_singlepass_342A/backgrounds'
angles = [280,297,312,322]
angle0 = 280

#raw data plots
fig, axs = plt.subplots(1, 2, figsize=(14, 8))

axs[0].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[0].set_ylabel("Single Beam",fontsize=12)
theta_variation_title =sample_name + ' Raw data'
axs[0].set_title(theta_variation_title)

fit_plot_title = sample_name+ " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax) + r" ${cm}^{-1}$"
axs[1].set_title(fit_plot_title)
axs[1].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[1].set_ylabel('Transmission Ratio',fontsize=12)

tm_bg_file = os.path.join(bg_dir, 'background_P0deg' + '.CSV')
te_bg_file = os.path.join(bg_dir, 'background_P90deg' + '.CSV')

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

#absorption fit plot
fig_fits, axs_fits = plt.subplots(figsize=(12, 8))
axs_fits.set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs_fits.set_ylabel(r"$\alpha_{ISB} \times L_{path}=-\ln (\frac{I_{out,TM}}{I_{out,TE}}) + \ln(\frac{I_{bg,TM}}{I_{bg,TE}})$",fontsize=12)
fit_plot_title = sample_name + " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax)
axs_fits.set_title(fit_plot_title)

fig_trans, axs_trans = plt.subplots(figsize=(12, 8))
axs_trans.set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs_trans.set_ylabel(r"$\frac{I_{out,TM(TE)}}{I_{in,TM(TE)}}$",fontsize=12)
fit_trans_title = sample_name + " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax)
axs_trans.set_title(fit_trans_title)

blues = ['darkblue','mediumblue','blue','cornflowerblue']
reds = ['mediumvioletred','deeppink','hotpink','palevioletred']

for i in range(0,len(angles)):
    angle = angles[i]
    angle_filename = 'theta_stage_'+ str(angle) + '_deg_'
    tm_file = os.path.join(base_dir, angle_filename+'P0deg' + '.CSV')
    te_file = os.path.join(base_dir, angle_filename+'P90deg' + '.CSV')
    _, tm_wavenum, tm_single_beam, _ = load_data(tm_file)
    _, te_wavenum, te_single_beam, _ = load_data(te_file)
    angle = angle-angle0
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
                            color=blues[i])

    #calculate the absorption coefficient times path length

    offset = np.log(bg_meas.TM_masked/bg_meas.TE_masked)

    alpha_ISB= -np.log(angle_meas.TM_masked/angle_meas.TE_masked)+offset
    axs_fits.plot(angle_meas.TE_wavenum_masked, alpha_ISB, label= str(angle) + '$\degree$',
                            color=blues[i])

    if angle == 0:
        axs_trans.plot(angle_meas.TM_wavenum_masked, angle_meas.TM_masked/bg_meas.TM_masked, label='TM '+ str(angle) + '$\degree$',
                            color=blues[i])
        axs_trans.plot(angle_meas.TE_wavenum_masked, angle_meas.TE_masked/bg_meas.TE_masked, label= 'TE '+str(angle) + '$\degree$',
                            color=reds[i])

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

#
plt.figure(fig_fits)
axs_fits.legend()
axs_fits.legend(prop={"size":14})
plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'alpha_ISBs' + '.svg')
plt.savefig(save_title)

plt.figure(fig_trans)
axs_trans.legend()
axs_trans.legend(prop={"size":14})
plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'transmissions' + '.svg')
plt.savefig(save_title)

plt.figure(fig)
plt.show()
plt.figure(fig_trans)
plt.show()
plt.figure(fig_fits)
plt.show()

