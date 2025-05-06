import os
import matplotlib.pyplot as plt
import numpy as np
from FTIR_analysis_helpers import load_data
from FTIR_analysis_helpers import SinglePassMeas
from FTIR_analysis_helpers import build_SP
from FTIR_analysis_helpers import calculate_Fresnels
from matplotlib.ticker import MaxNLocator

sample_name = 'P4-12-19-2A-beta1'
numax = 4000
numin = 650
ap_meas = 40
gain_meas = 4
ap_bg = 10
gain_bg = 2
Fout = 0.78 #fresnel coefficient in
Fin = Fout
nuISBs_check = [992, 1342, 1517]

Lsample = 10 #mm
tsamp = np.mean([0.449,0.443,0.445,0.438]) #mm
Nbounces = Lsample/tsamp
Lpath = 14

# settings_suffix_meas = 'ap-'+str(ap_meas)+'-gain-'+str(gain_meas)
# settings_suffix_bg = 'ap-'+str(ap_bg)+'-gain-'+str(gain_bg)

#adjust with well thicknesses based on Lodo runsheet

base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250505_P4-12-19-2A-beta1'
n_air = 1.0
n_InP = 2.7132  # InP at lambda = 8 um
n_GaAs = 3.28
angles = [0,30,50]

# angle0 = 280

#raw data plots
fig, axs = plt.subplots(1, 3, figsize=(14, 8))

axs[0].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[0].set_ylabel(r"$\frac{Single Beam}{T_{12} \times T_{23}}$",fontsize=12)
theta_variation_title =sample_name + ' Raw data, Fresnel corrected'
axs[0].set_title(theta_variation_title)

fit_plot_title = sample_name+ " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax) + r" ${cm}^{-1}$"
axs[1].set_title(fit_plot_title)
axs[1].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[1].set_ylabel('Transmission Ratio',fontsize=12)

fit_plot_title = "angle resolved transmissions fresnel corrected"
axs[2].set_title(fit_plot_title)
axs[2].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[2].set_ylabel('Transmission Ratio',fontsize=12)

axs[0].xaxis.set_major_locator(MaxNLocator(integer=True))
axs[1].xaxis.set_major_locator(MaxNLocator(integer=True))
axs[2].xaxis.set_major_locator(MaxNLocator(integer=True))

tm_bg_file = os.path.join(base_dir, 'bg' +'_P0deg' + '.CSV')
te_bg_file = os.path.join(base_dir, 'bg' + '_P90deg' + '.CSV')

bg_meas = build_SP(te_bg_file, tm_bg_file, 'bg', [0], fresnel=False, nuextrema=[numin, numax])

axs[0].plot(bg_meas.TE_wavenum, bg_meas.TE_single_beam, label='TE bg ',
            color='yellowgreen')
axs[0].plot(bg_meas.TM_wavenum, bg_meas.TM_single_beam, label='TM bg',
            color='darkgoldenrod')

axs[1].plot(bg_meas.TM_wavenum_masked, bg_meas.TM_masked/bg_meas.TE_masked,
            label='TM/TE no sample', color='y')

#absorption fit plot
fig_fits, axs_fits = plt.subplots(figsize=(12, 8))
axs_fits.set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs_fits.set_ylabel(r"$\alpha_{ISB} \times L_{path}=-\ln (\frac{I_{out,TM}/F}{I_{out,TE}/F}) + \ln(\frac{I_{bg,TM}}{I_{bg,TE}})$",fontsize=12)
fit_plot_title = sample_name + " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax)
axs_fits.set_title(fit_plot_title)
axs_fits.grid()
axs_fits.xaxis.set_major_locator(MaxNLocator(integer=True))


#absorption per pass plot
fig_abs, axs_abs = plt.subplots(1,2,figsize=(12, 8))
axs_abs[0].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs_abs[0].set_ylabel(r"$\alpha = \frac{-\ln(\frac{I_{out}}{I_{in}} \frac{1}{F_{in} \times F_{out}})}{t_{samp}} $",fontsize=12)
fit_abs_plot_title = sample_name + " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax) + r", $t_{samp} = $" + str(tsamp) + " mm"
axs_abs[0].set_title(fit_abs_plot_title)

axs_abs[1].set_xlabel(r"path length [mm]",fontsize=12)
axs_abs[1].set_ylabel(r"Transmission fraction per sample path length ",fontsize=12)

alpha_nu0_TE = 0 # per mm
alpha_nu0_TM = 0 #place holders to be filled in later

blues = ['darkblue','mediumblue','blue','cornflowerblue']
reds = ['fuchsia','mediumvioletred','crimson','lightpink']

for i in range(0,len(angles)):
    angle = angles[i]

    angle_filename = 'rot'+ str(angle)+'deg_'

    tm_file = os.path.join(base_dir, angle_filename+'P0deg' + '.CSV')
    te_file = os.path.join(base_dir, angle_filename+'P90deg' + '.CSV')

    angle_meas = build_SP(te_file,tm_file,sample_name,[angle],fresnel=True,n1=n_air,n2=n_GaAs,n3=n_air,nuextrema=[numin,numax])

    #add the raw data plot
    axs[0].plot(angle_meas.TE_wavenum, angle_meas.TE_single_beam, label='TE '+str(angle) + '$\degree$',
                color=blues[i])
    axs[0].plot(angle_meas.TM_wavenum, angle_meas.TM_single_beam, label='TM '+str(angle) + '$\degree$',
                color=reds[i])

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

    if angle==0:

        alpha_samp_TE = -np.log(angle_meas.TE_masked/bg_meas.TE_masked)/tsamp
        alpha_samp_TM = -np.log(angle_meas.TM_masked / bg_meas.TM_masked)/tsamp
        # nuISBs_check = [0.95*nu_ISB,nu_ISB,1.05*nu_ISB]

        axs_abs[0].plot(angle_meas.TE_wavenum_masked,alpha_samp_TE,label='TE' + str(angle) + '$\degree$',color=reds[i])
        axs_abs[0].plot(angle_meas.TM_wavenum_masked, alpha_samp_TM,label='TM' + str(angle) + '$\degree$',color=blues[i])
        for nuISB_ind in range(0,len(nuISBs_check)):
            nuISB_check = nuISBs_check[nuISB_ind]
            idx_closest = np.abs(angle_meas.TE_wavenum_masked - nuISB_check).argmin()

            alpha_nu0_TE = alpha_samp_TE[idx_closest] #per mm
            alpha_nu0_TM = alpha_samp_TM[idx_closest]
            pathlens = np.linspace(0.0,Lpath,num=50)#per mm
            axs_abs[1].plot(pathlens, np.exp(-alpha_nu0_TM*pathlens),label=r"$TM, \nu=$" + str(nuISB_check) + r"${cm}^{-1}$",color=reds[nuISB_ind])
            axs_abs[1].plot(pathlens, np.exp(-alpha_nu0_TE * pathlens), label=r"$TE, \nu=$" + str(nuISB_check) + r"${cm}^{-1}$",color=blues[nuISB_ind])

    #transmission per length at 1120 cm^-1


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
#                    color='green')

plt.figure(fig_fits)
axs_fits.legend()
axs_fits.legend(prop={"size":14})
plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'alpha_ISBs' + '.svg')
plt.savefig(save_title)

plt.figure(fig_abs)
# fit_abs_plot_title = r"$\alpha(\nu = $" + str(nu_ISB) + r")"
# axs_abs[1].set_title(fit_abs_plot_title)

axs_abs[0].grid(True)
axs_abs[1].grid(True)
axs_abs[0].legend()
axs_abs[0].legend(prop={"size":14})
axs_abs[1].legend()
axs_abs[1].legend(prop={"size":14})
plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'transmissions' + '.svg')
plt.savefig(save_title)

plt.figure(fig)
plt.show()
plt.figure(fig_fits)
plt.show()
plt.figure(fig_abs)
plt.show()