import os
import matplotlib.pyplot as plt
import numpy as np
from FTIR_analysis_helpers import load_data
from FTIR_analysis_helpers import SinglePassMeas

sample_name = 'P530-C'
numax = 3100
numin = 800
ap = 20
gain = 2
Fout = 0.78 #fresnel coefficient in
Fin = Fout
tsamp = 0.4 #mm
nu_ISB = 1120

settings_suffix = 'ap-'+str(ap)+'-gain-'+str(gain)

#adjust with well thicknesses based on Lodo runsheet

base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250204_P530-C-SP'
angles = [0]

# angle0 = 280

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

tm_bg_file = os.path.join(base_dir, 'bg-P0deg-' + settings_suffix + '.CSV')
te_bg_file = os.path.join(base_dir, 'bg-P90deg-' + settings_suffix + '.CSV')

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

#absorption per pass plot
fig_abs, axs_abs = plt.subplots(1,2,figsize=(12, 8))
axs_abs.set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs_abs.set_ylabel(r"$\alpha = \frac{-\ln(\frac{I_{out}}{I_{in}} \frac{1}{F_{in} \times F_{out}})}{t_{samp}} $",fontsize=12)
fit_abs_plot_title = sample_name + " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax)
axs_abs.set_title(fit_abs_plot_title)

blues = ['darkblue','mediumblue','blue','cornflowerblue']
reds = ['mediumvioletred','deeppink','hotpink','pink']

for i in range(0,len(angles)):
    angle = angles[i]
    angle_filename = sample_name + '-SP-rot-'+ str(angle)+'deg'

    tm_file = os.path.join(base_dir, angle_filename+'-P0deg-' + settings_suffix + '.CSV')
    te_file = os.path.join(base_dir, angle_filename+'-P90deg-' +settings_suffix + '.CSV')
    _, tm_wavenum, tm_single_beam, _ = load_data(tm_file)
    _, te_wavenum, te_single_beam, _ = load_data(te_file)
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

    #calculate the absorption coefficient times path length

    offset = np.log(bg_meas.TM_masked/bg_meas.TE_masked)

    alpha_ISB= -np.log(angle_meas.TM_masked/angle_meas.TE_masked)+offset
    axs_fits.plot(angle_meas.TE_wavenum_masked, alpha_ISB, label= str(angle) + '$\degree$',
                            color=reds[i])

    axs_abs[0].plot(angle_meas.TE_wavenum_masked,-np.log((angle_meas.TE_masked/bg_meas.TE_masked)/(Fin*Fout))/tsamp,label='TE' + str(angle) + '$\degree$')
    axs_abs[0].plot(angle_meas.TM_wavenum_masked, -np.log((angle_meas.TM_masked / bg_meas.TM_masked)/(Fin*Fout))/tsamp,label='TM' + str(angle) + '$\degree$')

    #transmission per length at 1120 cm^-1


plt.figure(fig)
axs[0].legend()
axs[0].legend(prop={"size":14})
axs[1].legend()
axs[1].legend(prop={"size":14})
plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'raw scans and ratios' + '.svg')
plt.savefig(save_title)
#                    color='green')

plt.figure(fig_fits)
axs_fits.legend()
axs_fits.legend(prop={"size":14})
plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'alpha_ISBs' + '.svg')
plt.savefig(save_title)

plt.figure(fig_abs)
axs_abs.legend()
axs_abs.legend(prop={"size":14})
plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'transmissions' + '.svg')
plt.savefig(save_title)

plt.figure(fig)
plt.show()
plt.figure(fig_fits)
plt.show()
plt.figure(fig_abs)
plt.show()