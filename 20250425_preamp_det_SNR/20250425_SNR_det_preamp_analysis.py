import os
# import scipy
import matplotlib.pyplot as plt
import numpy as np
from FTIR_analysis_helpers import MultipassMeas
from FTIR_analysis_helpers import load_data
from matplotlib.ticker import MaxNLocator
from FTIR_analysis_helpers import fitLorentzPlot

ap = 15
gain = 2
att = 'mod'
preamps = ['04-07-101546','19-05-102110']

dets = ['J15D22-M204-S01M-60','J15D16-S']

path = ['P4-12-19-2A-alpha','none']


numax = 3900
numin = 600

#dets to load
#preamps to load

#adjust with well thicknesses based on Lodo runsheet

base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/preamp_SNR_test/none_in_path'

fig, axs = plt.subplots(1, 3, figsize=(14, 8))

overall_title ='SNR comparison'
fig.suptitle(overall_title)

axs[0].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[0].set_ylabel("Single Beam",fontsize=12)
theta_variation_title ='Raw data ' + dets[0] + ' att = none, ' + dets[1] + ' att = moderate'
axs[0].set_title(theta_variation_title)

fit_plot_title = "mask " + str(numin) + r"$ < \nu < $" + str(numax) + r" ${cm}^{-1}$"
axs[1].set_title(fit_plot_title)
axs[1].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[1].set_ylabel('Transmission with ' + preamps[1] + '/ Transmission with ' + preamps[0] ,fontsize=12)

fit_plot_title = "noise fits"
axs[2].set_title(fit_plot_title)
axs[2].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[2].set_ylabel('Single Beam',fontsize=12)

meas_dict = dict()

colors = ['red','dodgerblue','aquamarine','fuchsia']

color_count = 0

for preamp in preamps:
    for det in dets:
        #load file

        #create struct
        #add to plot
        meas_name = 'preamp ' + preamp + ' det ' + det
        meas = MultipassMeas(samp=meas_name)
        filename = det + '-det-SN-' + preamp + '-preamp'
        fileloaded = os.path.join(base_dir, filename + '.CSV')

        _, wavenum, single_beam, _ = load_data(fileloaded)

        meas.TE_single_beam = single_beam
        meas.TE_wavenum = wavenum

        mask_samp = (meas.TE_wavenum > numin) & (meas.TE_wavenum < numax)
        meas.TE_masked = meas.TE_single_beam[mask_samp]

        meas.TE_wavenum_masked = meas.TE_wavenum[mask_samp]

        axs[0].plot(meas.TE_wavenum, meas.TE_single_beam, color=colors[color_count], label=meas.name)

        meas_dict[meas.name] = meas

        #add the polynomial fit
        mask_std = (meas.TE_wavenum > 2690) & (meas.TE_wavenum < 2720)

        slope, intercept = np.polyfit(meas.TE_wavenum[mask_std], meas.TE_single_beam[mask_std], 1)
        yfit = slope * meas.TE_wavenum[mask_std] + intercept
        std = np.std(meas.TE_single_beam[mask_std]-yfit)

        # axs[2].scatter(meas.TE_wavenum[mask_std], meas.TE_single_beam[mask_std], label=meas.name + 'raw',s=1)
        axs[2].scatter(meas.TE_wavenum[mask_std], meas.TE_single_beam[mask_std],color=colors[color_count],s=1)

        axs[2].plot(meas.TE_wavenum[mask_std], yfit, color=colors[color_count],label="std = %0.4f" % std)

        color_count = color_count+1
#try using backgrounds with no sample in the path
for det in dets:
    axs[1].plot(meas_dict['preamp ' + preamps[0] + ' det ' + det].TE_wavenum_masked, meas_dict['preamp ' + preamps[1] + ' det ' + det].TE_masked/meas_dict['preamp ' + preamps[0] + ' det ' + det].TE_masked, label=det)


plt.figure(fig)
axs[0].legend()
axs[0].legend(prop={"size":14})
axs[1].legend()
axs[1].legend(prop={"size":14})
axs[2].legend()
axs[2].legend(prop={"size":14})
plt.tight_layout()
save_title = os.path.join(base_dir, overall_title + 'raw scans and ratios' + '.svg')
plt.savefig(save_title)

plt.figure(fig)
plt.show()