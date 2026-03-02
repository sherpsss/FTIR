import os
# import scipy
import matplotlib.pyplot as plt
import numpy as np
from FTIR_analysis_helpers import MultipassMeas
from FTIR_analysis_helpers import load_data
from FTIR_analysis_helpers import build_MP
from matplotlib.ticker import MaxNLocator
from FTIR_analysis_helpers import fitLorentzPlot

sample_name = 'P12-9-25-1-GaAs-intrinsic-only-MP'

numax = 3100
numin = 820

base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/P12-9-25-1-intrinsic-GaAs-only-MP'
bg_dir = base_dir

tm_file = os.path.join(base_dir, sample_name+ '-P0deg' + '.CSV')
te_file = os.path.join(base_dir, sample_name+ '-P90deg' + '.CSV')

samp_meas = build_MP(te_file,tm_file,sample_name,nuextrema=[numin,numax])

tm_bg_file = os.path.join(bg_dir,'P0deg_' + 'bg' + '.CSV')
te_bg_file = os.path.join(bg_dir, 'P90deg_' + 'bg' + '.CSV')

bg_meas = build_MP(te_bg_file,tm_bg_file,'no-samp',nuextrema=[numin,numax])

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
axs[2].set_ylabel(samp_meas.name + ' / ' + bg_meas.name + ' transmission ratios',fontsize=12)

axs[0].grid()
axs[0].xaxis.set_major_locator(MaxNLocator(integer=True))
axs[1].grid()
axs[1].xaxis.set_major_locator(MaxNLocator(integer=True))
axs[2].grid()
axs[2].xaxis.set_major_locator(MaxNLocator(integer=True))

axs[0].plot(samp_meas.TE_wavenum, samp_meas.TE_single_beam, label= 'TE' + samp_meas.name,
                        color='blue')
axs[0].plot(samp_meas.TM_wavenum , samp_meas.TM_single_beam , label= 'TM' + samp_meas.name,
                        color='red')

axs[0].plot(bg_meas.TE_wavenum, bg_meas.TE_single_beam, label= 'TE ' + bg_meas.name,
                        color='c')
axs[0].plot(bg_meas.TM_wavenum , bg_meas.TM_single_beam , label='TM ' + bg_meas.name,
                        color='m')


axs[1].plot(samp_meas.TM_wavenum_masked, samp_meas.TM_masked / samp_meas.TE_masked,
            label='TM/TE ' + sample_name)
axs[1].plot(bg_meas.TM_wavenum_masked, bg_meas.TM_masked / bg_meas.TE_masked,
            label='TM/TE ' + 'bg')

axs[2].plot(samp_meas.TM_wavenum_masked, samp_meas.TM_masked / bg_meas.TM_masked,
            label='TM ' + samp_meas.name + '/ TM ' + bg_meas.name)
axs[2].plot(samp_meas.TE_wavenum_masked, samp_meas.TE_masked / bg_meas.TE_masked,
            label='TE ' + samp_meas.name + ' / TE ' + bg_meas.name)

plt.figure(fig)
axs[0].legend()
axs[0].legend(prop={"size":10})
axs[1].legend()
axs[1].legend(prop={"size":10})
axs[2].legend()
axs[2].legend(prop={"size":10})
plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'raw scans and ratios' + '.svg')
plt.savefig(save_title)

fig_fits, axs_fits = plt.subplots(figsize=(10, 8))
axs_fits.set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs_fits.set_ylabel(r"$\alpha_{ISB} \times L{path} $",fontsize=12)
fit_plot_title = sample_name + " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax)
axs_fits.set_title(fit_plot_title)
axs_fits.grid()
axs_fits.xaxis.set_major_locator(MaxNLocator(integer=True))
#
offset = np.log(bg_meas.TM_masked/bg_meas.TE_masked)

alpha_ISB= -np.log(samp_meas.TM_masked/samp_meas.TE_masked)+offset

axs_fits.plot(
    samp_meas.TM_wavenum_masked,
    alpha_ISB,
    color='green',
    label=rf"$-\ln \left(\frac{{I_{{{samp_meas.name},TM}}}}{{I_{{{samp_meas.name},TE}}}}\right)"
          rf"+ \ln \left(\frac{{I_{{{bg_meas.name},TM}}}}{{I_{{{bg_meas.name},TE}}}}\right)$"
)

# numins = [1060]
# numaxs = [1670]
# kappa_nu_guesses = [65]
# for i in range(0,len(numins)):
#     nu_range = [numins[i],numaxs[i]]
#     kappa_nu_guess = kappa_nu_guesses[i]
#     fitLorentzParams = fitLorentzPlot(nu_range, kappa_nu_guess,samp_meas.TE_wavenum_masked, alpha_ISB, axs_fits)

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