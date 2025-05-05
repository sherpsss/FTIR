import os
# import scipy
import matplotlib.pyplot as plt
import numpy as np
from FTIR_analysis_helpers import MultipassMeas
from FTIR_analysis_helpers import load_data
from FTIR_analysis_helpers import build_MP
from matplotlib.ticker import MaxNLocator
from FTIR_analysis_helpers import fitLorentzPlot



sample_name = '342B'
ap_bg = 15
ap_samp = 87
gain_bg = 2
gain_samp = 4
att_bg = 'mod'
preamp = '04-07-101546'

# settings_suffix_bg = 'ap_' + str(ap_bg) + '_gain_' + str(gain_bg)
settings_suffix_bg = 'att_' + att_bg

settings_suffix_samp = 'ap_' + str(ap_samp) + '_gain_' + str(gain_samp)

numax = 3400
numin = 660


#adjust with well thicknesses based on Lodo runsheet


base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250422_342B'
bg_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250421_P4-12-19-2A-alpha/PA101-J15D16-preamp-matched'

# tm_file = os.path.join(base_dir, settings_suffix_samp + '-P0deg' + '.CSV')
# te_file = os.path.join(base_dir, settings_suffix_samp + '-P90deg' + '.CSV')
tm_file = os.path.join(base_dir, 'P0deg_' + settings_suffix_samp + '.CSV')
te_file = os.path.join(base_dir, 'P90deg_' + settings_suffix_samp + '.CSV')

samp_meas = build_MP(te_file,tm_file,sample_name,nuextrema=[numin,numax])

tm_bg_file = os.path.join(bg_dir,'bg_'+settings_suffix_bg+ '_P0deg' + '.CSV')
te_bg_file = os.path.join(bg_dir, 'bg_'+settings_suffix_bg+'_P90deg' + '.CSV')

bg_meas = build_MP(te_bg_file,tm_bg_file,'no_samp',nuextrema=[numin,numax])

# bg_name = 'no_samp'
# bg_meas = MultipassMeas(samp=bg_name)


# _, tm_bg_wavenum, tm_bg_single_beam, _ = load_data(tm_bg_file)
# _, te_bg_wavenum, te_bg_single_beam, _ = load_data(te_bg_file)
#
# bg_meas.TE_single_beam=te_bg_single_beam
# bg_meas.TM_single_beam=tm_bg_single_beam
# bg_meas.TE_wavenum=te_bg_wavenum
# bg_meas.TM_wavenum=tm_bg_wavenum

fig, axs = plt.subplots(1, 3, figsize=(14, 8))

overall_title =sample_name
fig.suptitle(overall_title)

axs[0].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[0].set_ylabel("Single Beam",fontsize=12)
theta_variation_title =' Raw data'
axs[0].set_title(theta_variation_title)

fit_plot_title = " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax) + r" ${cm}^{-1}$"
axs[1].set_title(fit_plot_title)
axs[1].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[1].set_ylabel('Polarization Transmission Ratio',fontsize=12)

fit_plot_title = " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax) + r" ${cm}^{-1}$"
axs[2].set_title(fit_plot_title)
axs[2].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[2].set_ylabel('Samp over bg transmission ratios',fontsize=12)

axs[0].grid()
axs[0].xaxis.set_major_locator(MaxNLocator(integer=True))
axs[1].grid()
axs[1].xaxis.set_major_locator(MaxNLocator(integer=True))
axs[2].grid()
axs[2].xaxis.set_major_locator(MaxNLocator(integer=True))

axs[0].plot(samp_meas.TE_wavenum, samp_meas.TE_single_beam, label= 'TE',
                        color='blue')
axs[0].plot(samp_meas.TM_wavenum , samp_meas.TM_single_beam , label= 'TM',
                        color='red')

axs[0].plot(bg_meas.TE_wavenum, bg_meas.TE_single_beam, label= 'TE bg ' + 'no_samp',
                        color='c')
axs[0].plot(bg_meas.TM_wavenum , bg_meas.TM_single_beam , label='TM bg' + 'no samp',
                        color='m')

# do the masking
# mask_samp = (bg_meas.TE_wavenum > numin) & (bg_meas.TE_wavenum < numax)

# bg_meas.TE_masked = bg_meas.TE_single_beam[mask_samp]
# bg_meas.TM_masked = bg_meas.TM_single_beam[mask_samp]
# bg_meas.TM_wavenum_masked = bg_meas.TM_wavenum[mask_samp]
# bg_meas.TE_wavenum_masked = bg_meas.TE_wavenum[mask_samp]

#plot the masked signals
#
# axs[0].plot(samp_meas.TE_wavenum_masked, samp_meas.TE_masked, label= 'TE' + ' masked',
#                         color='cornflowerblue')
# axs[0].plot(samp_meas.TM_wavenum_masked , samp_meas.TM_masked , label= 'TM' + ' masked',
#                         color='tomato')
#
# axs[0].plot(bg_meas.TE_wavenum_masked, bg_meas.TE_masked, label= 'TE bg ' + bg_name + ' masked',
#                         color='cadetblue')
# axs[0].plot(bg_meas.TM_wavenum_masked , bg_meas.TM_masked , label='TM bg ' + bg_name + ' masked',
#                         color='indigo')


axs[1].plot(samp_meas.TM_wavenum_masked, samp_meas.TM_masked / samp_meas.TE_masked,
            label='TM/TE ' + sample_name)
axs[1].plot(bg_meas.TM_wavenum_masked, bg_meas.TM_masked / bg_meas.TE_masked,
            label='TM/TE ' + 'bg')

axs[2].plot(samp_meas.TM_wavenum_masked, samp_meas.TM_masked / bg_meas.TM_masked,
            label='TM samp/ TM bg')
axs[2].plot(samp_meas.TE_wavenum_masked, samp_meas.TE_masked / bg_meas.TE_masked,
            label='TE samp/ TE bg')

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

fig_fits, axs_fits = plt.subplots(figsize=(10, 8))
axs_fits.set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs_fits.set_ylabel(r"$-\ln (\frac{I_{out,TM}}{I_{out,TE}}) + \ln(\frac{I_{bg,TM}}{I_{bg,TE}})$",fontsize=12)
fit_plot_title = sample_name + " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax)
axs_fits.set_title(fit_plot_title)
axs_fits.grid()
axs_fits.xaxis.set_major_locator(MaxNLocator(integer=True))
#
offset = np.log(bg_meas.TM_masked/bg_meas.TE_masked)

alpha_ISB= -np.log(samp_meas.TM_masked/samp_meas.TE_masked)+offset
#
# # axs_fits.plot(samp_meas.TE_wavenum_masked, alpha_ISB, label= r"$\alpha_{ISB}$ per well",
#                         # color='green')
axs_fits.plot(samp_meas.TM_wavenum_masked, alpha_ISB, color='green')

plt.figure(fig_fits)
axs_fits.legend()
axs_fits.legend(prop={"size":14})

plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'alpha_ISB' + '.svg')
plt.savefig(save_title)

plt.figure(fig)
plt.show()
plt.figure(fig_fits)
plt.show()