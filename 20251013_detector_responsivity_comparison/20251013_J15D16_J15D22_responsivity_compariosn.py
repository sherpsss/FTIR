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

base_dir_J15D16 = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20251013_J15D16_det'
base_dir_J15D22 = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20251014_J15D22_in_path'
no_samp_path = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20241114_unpat_P551/no samp in path'
n_air = 1.0
n_InP = 2.7132  # InP at lambda = 8 um
n_GaAs = 3.28
detectors = ['J15D22','J15D16']
paths = [base_dir_J15D22,base_dir_J15D16]

# angle0 = 280

#raw data plots

#plot all the raw data
fig, axs = plt.subplots(1, 2, figsize=(18, 8))

axs[0].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[0].set_ylabel(r"Single Beam",fontsize=12)
theta_variation_title =sample_name + ' Raw data'
axs[0].set_title(theta_variation_title)

# axs[1].set_title(fit_plot_title)
axs[1].set_title('normalized')
axs[1].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[1].set_ylabel('single beam/max(single beam)',fontsize=12)

axs[1].xaxis.set_major_locator(MaxNLocator(integer=True))
axs[1].xaxis.set_major_locator(MaxNLocator(integer=True))

blues = ['darkblue','mediumblue','blue','cornflowerblue']
reds = ['fuchsia','mediumvioletred','crimson','lightpink']
greens = ['seagreen','mediumseagreen']
# _, wavenum_np, sb_np, _ = load_data(nopol_file)
# sb_np_norm = sb_np / np.nanmax(sb_np)
# axs.plot(wavenum_np, sb_np_norm, '--', color='gray',
#                  label=f'No polarizer, det = {det_tag}')
# axs.plot(angle_meas.TE_wavenum, angle_meas.TE_single_beam_raw, label='TE, det = ' + det_tag,
#          color=blues[i])

for i in range(0,len(detectors)):
    det_tag = detectors[i]
    tm_file = None
    te_file = None
    nopol_file = None

    for fname in os.listdir(paths[i]):
        f_lower = fname.lower()  # make it case-insensitive
        if all(tag in f_lower for tag in ['p0deg', det_tag.lower(), '.csv']):
            tm_file = os.path.join(paths[i], fname)
        elif all(tag in f_lower for tag in ['p90deg', det_tag.lower(), '.csv']):
            te_file = os.path.join(paths[i], fname)
        elif all(tag in f_lower for tag in ['nopolarizer', det_tag.lower(), '.csv']):
            nopol_file = os.path.join(paths[i], fname)
    if det_tag == 'J15D16':
        for fname in os.listdir(no_samp_path):
            f_lower = fname.lower()
            if all(tag in f_lower for tag in ['nopolarizer','.csv']):
                nopol_file = os.path.join(no_samp_path,fname)

    if tm_file is None:
        raise FileNotFoundError(f"No TM file found in {paths[i]} with {det_tag} and P0deg")
    if te_file is None:
        raise FileNotFoundError(f"No TE file found in {paths[i]} with {det_tag} and P90deg")
    if nopol_file is None:
        print(f"⚠️ No 'no polarizer' file found for {det_tag}")
    else:
        print("No polarizer file:", nopol_file)

    print("TM file:", tm_file)
    print("TE file:", te_file)


    angle_meas = build_SP(te_file,tm_file,sample_name,[det_tag],fresnel=False,n1=n_air,n2=n_GaAs,n3=n_air,nuextrema=[numin,numax])

    #add the raw data plot
    axs[0].plot(angle_meas.TE_wavenum, angle_meas.TE_single_beam_raw, label='TE, det = '+det_tag,
                color=blues[i])
    axs[0].plot(angle_meas.TM_wavenum, angle_meas.TM_single_beam_raw, label='TM, ap =  '+det_tag,
                color=reds[i])

    axs[1].plot(angle_meas.TE_wavenum, angle_meas.TE_single_beam_raw/np.max(angle_meas.TE_single_beam_raw), label= 'TE, det = '+det_tag,
                            color=blues[i])
    axs[1].plot(angle_meas.TM_wavenum,angle_meas.TM_single_beam_raw / np.max(angle_meas.TM_single_beam_raw),
                     label='TM, det = ' + det_tag,color=reds[i])

    if nopol_file is not None:
        _, wavenum_np, singlebeam_np, _ = load_data(nopol_file)
        sb_np_norm = singlebeam_np / np.max(singlebeam_np)
        axs[0].plot(wavenum_np, singlebeam_np, label=f'No polarizer, det =' + det_tag, color=greens[i])
        axs[1].plot(wavenum_np,sb_np_norm,label='no polarizer, det = ' + det_tag, color=greens[i])


plt.figure(fig)
axs[0].grid()
axs[1].grid()
axs[0].legend()
axs[0].legend(prop={"size":14})
plt.tight_layout()
plt.tight_layout()
save_title = os.path.join(paths[1], sample_name + 'raw scans and ratios' + '.svg')
plt.savefig(save_title)
#                    color='green'


plt.figure(fig)
plt.show()