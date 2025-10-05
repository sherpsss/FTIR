import os
import matplotlib.pyplot as plt
import numpy as np
from FTIR_analysis_helpers import load_data
from FTIR_analysis_helpers import SinglePassMeas
from FTIR_analysis_helpers import build_SP
from FTIR_analysis_helpers import calculate_Fresnels
from matplotlib.ticker import MaxNLocator

sample_name = 'no_samp'
numax = 3300
numin = 700
gain_ap_sweep = 1

# settings_suffix_meas = 'ap-'+str(ap_meas)+'-gain-'+str(gain_meas)
# settings_suffix_bg = 'ap-'+str(ap_bg)+'-gain-'+str(gain_bg)

#adjust with well thicknesses based on Lodo runsheet

base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250922_FTIR1_ap_gain_dispersion/ap_sweep'
n_air = 1.0
n_InP = 2.7132  # InP at lambda = 8 um
n_GaAs = 3.28
aps = [25,50,100]

# angle0 = 280

#raw data plots
fig, axs = plt.subplots(1, 1, figsize=(18, 8))

fig_rats, axs_rats = plt.subplots(1, 2, figsize=(18, 8))

axs.set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs.set_ylabel(r"Single Beam",fontsize=12)
theta_variation_title =sample_name + ' Raw data'
axs.set_title(theta_variation_title)

fit_plot_title = sample_name+ " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax) + r" ${cm}^{-1}$"
axs_rats[0].set_title(fit_plot_title)
axs_rats[0].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs_rats[0].set_ylabel('Transmission Ratio',fontsize=12)

axs_rats[0].xaxis.set_major_locator(MaxNLocator(integer=True))
axs_rats[0].xaxis.set_major_locator(MaxNLocator(integer=True))
axs.xaxis.set_major_locator(MaxNLocator(integer=True))
axs.xaxis.set_major_locator(MaxNLocator(integer=True))

blues = ['darkblue','mediumblue','blue','cornflowerblue']
reds = ['fuchsia','mediumvioletred','crimson','lightpink']

#create smallest aperture to normalize against
angle_filename = 'ap' + str(aps[0]) + '_gain' + str(gain_ap_sweep)

tm_file = os.path.join(base_dir, 'P0deg_' + angle_filename + '.CSV')
te_file = os.path.join(base_dir, 'P90deg_' + angle_filename + '.CSV')

min_ap = build_SP(te_file, tm_file, sample_name, [aps[0]], fresnel=False, n1=n_air, n2=n_GaAs, n3=n_air,
                      nuextrema=[numin, numax])

for i in range(0,len(aps)):
    ap = aps[i]

    angle_filename = 'ap'+ str(ap)+'_gain' + str(gain_ap_sweep)

    tm_file = os.path.join(base_dir,'P0deg_' + angle_filename + '.CSV')
    te_file = os.path.join(base_dir, 'P90deg_'+ angle_filename + '.CSV')

    angle_meas = build_SP(te_file,tm_file,sample_name,[ap],fresnel=False,n1=n_air,n2=n_GaAs,n3=n_air,nuextrema=[numin,numax])

    #add the raw data plot
    axs.plot(angle_meas.TE_wavenum, angle_meas.TE_single_beam_raw, label='TE, ap = '+str(ap),
                color=blues[i])
    axs.plot(angle_meas.TM_wavenum, angle_meas.TM_single_beam_raw, label='TM, ap =  '+str(ap),
                color=reds[i])

    axs_rats[0].plot(angle_meas.TM_wavenum_masked, angle_meas.TM_masked/angle_meas.TE_masked, label= 'TM/TE, ap = '+str(ap),
                            color=reds[i])

    if ap is not aps[0]:

        axs_rats[1].plot(angle_meas.TM_wavenum_masked, angle_meas.TM_masked/min_ap.TM_masked, label= 'TM ap = '+str(ap) + '/TM ap ' + str(aps[0]),
                            color=reds[i])
        axs_rats[1].plot(angle_meas.TE_wavenum_masked, angle_meas.TE_masked/min_ap.TE_masked, label= 'TE ap = '+str(ap) + '/TE ap ' + str(aps[0]),
                            color=blues[i])

    #transmission per length at 1120 cm^-1


plt.figure(fig)
axs.legend()
axs.legend(prop={"size":14})
plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'raw scans and ratios' + '.svg')
plt.savefig(save_title)
#                    color='green')

plt.figure(fig_rats)
axs_rats[0].legend()
axs_rats[0].legend(prop={"size":14})
axs_rats[1].legend()
axs_rats[1].legend(prop={"size":14})
plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + ' ratios' + '.svg')
plt.savefig(save_title)


plt.figure(fig)
plt.show()
plt.figure(fig_rats)
plt.show()