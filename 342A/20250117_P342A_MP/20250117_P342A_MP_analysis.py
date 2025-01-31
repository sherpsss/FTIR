import os
import matplotlib.pyplot as plt
import numpy as np
from FTIR_analysis_helpers import MultipassMeas
from FTIR_analysis_helpers import load_data

sample_name = 'P342A-MP'
wafername = 'P342A'
aperture = 87
gain = 8

numax = 3100
numin = 1450

blues = ['cornflowerblue','mediumblue','blue']
reds = ['crimson','m','magenta']

#adjust with well thicknesses based on Lodo runsheet

samp_meas = MultipassMeas(samp=sample_name)

base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250117_342A_multipass'

tm_file = os.path.join(base_dir, '342A_better_align_P0deg' + '.CSV')
te_file = os.path.join(base_dir, '342A_better_alignment_P90deg' + '.CSV')

_, tm_wavenum, tm_single_beam, _ = load_data(tm_file)
_, te_wavenum, te_single_beam, _ = load_data(te_file)

samp_meas.TM_wavenum = tm_wavenum
samp_meas.TE_wavenum = te_wavenum
samp_meas.TM_single_beam = tm_single_beam
samp_meas.TE_single_beam = te_single_beam

bg_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250122_P530-A-SP'

tm_bg_file = os.path.join(bg_dir, 'no_samp_P0deg' + '.CSV')
te_bg_file = os.path.join(bg_dir, 'no_samp_P90deg' + '.CSV')

_, tm_bg_wavenum, tm_bg_single_beam, _ = load_data(tm_bg_file)
_, te_bg_wavenum, te_bg_single_beam, _ = load_data(te_bg_file)

bg_meas = MultipassMeas(samp='background')

bg_meas.TE_single_beam=te_bg_single_beam
bg_meas.TM_single_beam=tm_bg_single_beam
bg_meas.TE_wavenum=te_bg_wavenum
bg_meas.TM_wavenum=tm_bg_wavenum

fig, axs = plt.subplots(1, 2, figsize=(14, 8))

axs[0].plot(samp_meas.TE_wavenum, samp_meas.TE_single_beam, label= 'TE MP',color=blues[0])
axs[0].plot(samp_meas.TM_wavenum , samp_meas.TM_single_beam , label= 'TM MP', color=reds[0])
axs[0].plot(bg_meas.TE_wavenum, bg_meas.TE_single_beam, label='TE bg from different day',
            color='yellowgreen')
axs[0].plot(bg_meas.TM_wavenum, bg_meas.TM_single_beam, label='TM bg from different day',
            color='darkgoldenrod')

mask_samp = (samp_meas.TE_wavenum > numin) & (samp_meas.TE_wavenum < numax)

# do the masking
samp_meas.TE_masked = samp_meas.TE_single_beam[mask_samp]
samp_meas.TM_masked = samp_meas.TM_single_beam[mask_samp]
samp_meas.TM_wavenum_masked = samp_meas.TM_wavenum[mask_samp]
samp_meas.TE_wavenum_masked = samp_meas.TE_wavenum[mask_samp]

mask_bg = (bg_meas.TE_wavenum > numin) & (bg_meas.TE_wavenum < numax)

bg_meas.TE_masked = bg_meas.TE_single_beam[mask_bg]
bg_meas.TM_masked = bg_meas.TM_single_beam[mask_bg]
bg_meas.TM_wavenum_masked = bg_meas.TM_wavenum[mask_bg]
bg_meas.TE_wavenum_masked = bg_meas.TE_wavenum[mask_bg]

axs[1].plot(samp_meas.TE_wavenum_masked, samp_meas.TM_masked/samp_meas.TE_masked, label= 'TE/TM masked',color=blues[0])
axs[1].plot(bg_meas.TE_wavenum_masked , bg_meas.TM_masked/bg_meas.TE_masked , label= 'bg ratio',color='yellowgreen')


axs[0].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[0].set_ylabel("Single Beam",fontsize=12)
# theta_variation_title =sample_name + ' compared to ' + SP_name
theta_variation_title = wafername + ' multipass, ap =  ' + str(aperture) + ', gain = ' + str(gain)
axs[0].set_title(theta_variation_title)

axs[0].legend()
axs[0].legend(prop={"size":14})

fit_plot_title = sample_name + " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax)
axs[1].set_title(fit_plot_title)
axs[1].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[1].set_ylabel('ratio',fontsize=12)

axs[1].legend()
axs[1].legend(prop={"size":14})


plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'raw data and ratios' + '.svg')
plt.savefig(save_title)
# plt.show()

#try checking to see if other ISBs are visible at the shorter wavelengths

fig_fits, axs_fits = plt.subplots(figsize=(12, 8))
axs_fits.set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs_fits.set_ylabel(r"$\alpha_{ISB} \times L_{path}=-\ln (\frac{I_{out,TM}}{I_{out,TE}}) + \ln(\frac{I_{bg,TM}}{I_{bg,TE}})$",fontsize=12)
fit_plot_title = sample_name + " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax)
axs_fits.set_title(fit_plot_title)

offset = np.log(bg_meas.TM_masked / bg_meas.TE_masked)

alpha_ISB = -np.log(samp_meas.TM_masked / samp_meas.TE_masked) + offset
axs_fits.plot(samp_meas.TE_wavenum_masked, alpha_ISB, label='absorption ',
              color=reds[1])

plt.figure(fig_fits)
axs_fits.legend()
axs_fits.legend(prop={"size":14})
plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'alpha_ISBs' + '.svg')
plt.savefig(save_title)

plt.figure(fig)
plt.show()
plt.figure(fig_fits)
plt.show()