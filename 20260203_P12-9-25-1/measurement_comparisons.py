import os
# import scipy
import matplotlib.pyplot as plt
import numpy as np
from FTIR_analysis_helpers import MultipassMeas
from FTIR_analysis_helpers import load_data
from FTIR_analysis_helpers import build_MP
from matplotlib.ticker import MaxNLocator
from FTIR_analysis_helpers import fitLorentzPlot


sample_name = 'P12-9-25-1-MP'
ap_bg = 15
gain_bg = 2
att_bg = 'mod'
preamp = 'PA101'

# settings_suffix_bg = 'ap_' + str(ap_bg) + '_gain_' + str(gain_bg)
settings_suffix_bg = 'att_' + att_bg

settings_suffix_samp = preamp + '-preamp-matched'

numax = 3100
numin = 825

base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20260203_P12-9-25-1-MP'
second_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20260203_P12-9-25-1-MP/20260206'
P129252MP_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20260205_P12-9-25-2-MP'
P129251SP_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20260204_P12-9-25-1-SP'

te_bg_file = os.path.join(base_dir, 'P90deg-bg' + '.CSV')
tm_bg_file = os.path.join(base_dir, 'P0deg-bg' + '.CSV')
te_bg_file_2nd = os.path.join(second_dir, 'P90deg_bg' + '.CSV')
tm_bg_file_2nd = os.path.join(second_dir,'P0deg_bg' + '.CSV')

tm_bg_file_P129252_MP = os.path.join(P129252MP_dir,'P0deg_bg' + '.CSV')
te_bg_file_P129252_MP = os.path.join(P129252MP_dir,'P90deg_bg' + '.CSV')

tm_bg_file_P129251SP = os.path.join(P129251SP_dir, 'P0deg_bg' + '.CSV')
te_bg_file_P129251SP = os.path.join(P129251SP_dir, 'P90deg_bg' + '.CSV')

fig, axs = plt.subplots(figsize=(14, 8))

overall_title =sample_name
fig.suptitle(overall_title)

axs.set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs.set_ylabel("Single Beam",fontsize=12)
theta_variation_title =' Raw data'
axs.set_title(theta_variation_title)

axs.grid()
axs.xaxis.set_major_locator(MaxNLocator(integer=True))

files = [te_bg_file,tm_bg_file,te_bg_file_2nd,tm_bg_file_2nd,te_bg_file_P129252_MP,tm_bg_file_P129252_MP,tm_bg_file_P129251SP,te_bg_file_P129251SP]

for file in files:
    # load the data
    print(file)
    wavelength, wavenum, single_beam, file_label,parsed_date = load_data(file,return_date=True)
    print(parsed_date)
    print(single_beam)
    # plot the data with the value given by the filename title
    axs.plot(wavenum, single_beam, label=file_label + str(parsed_date))
plt.figure(fig)
axs.legend()
axs.legend(prop={"size":14})

plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'bg comparisons' + '.svg')
plt.savefig(save_title)

plt.figure(fig)
plt.show()