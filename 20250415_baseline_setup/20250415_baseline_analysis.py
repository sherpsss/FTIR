import os
import matplotlib.pyplot as plt
import numpy as np
from FTIR_analysis_helpers import load_data
from FTIR_analysis_helpers import SinglePassMeas
from matplotlib.ticker import MaxNLocator

meas_date = '20250415'
baseline_date = '20250204'
ap = 20
gain = 2
settings_suffix = 'ap-'+str(ap)+'-gain-'+str(gain)

base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250415_setup_baseline/preamp-2ndstage'
baseline_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250204_P530-C-SP'

# angle0 = 280

#raw data plots
fig, axs = plt.subplots(1, 2, figsize=(14, 8))

axs[0].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[0].set_ylabel("Single Beam",fontsize=12)
overall_title =meas_date + 'vs.' + baseline_date + settings_suffix
fig.suptitle(overall_title)
axs[0].set_title('Raw data')

fit_plot_title = 'change in signal'
axs[1].set_title(fit_plot_title)
axs[1].set_xlabel('Wavenumber (cm^-1)',fontsize=12)
axs[1].set_ylabel(meas_date + ' Single beam/ ' + baseline_date + ' Single beam',fontsize=12)

axs[0].grid()
axs[1].grid()
# axs[2].grid()

axs[0].xaxis.set_major_locator(MaxNLocator(integer=True))

tm_bg_file = os.path.join(baseline_dir, 'bg-P0deg-' + settings_suffix + '.CSV')
te_bg_file = os.path.join(baseline_dir, 'bg-P90deg-' + settings_suffix+ '.CSV')

_, tm_bg_wavenum, tm_bg_single_beam, _ = load_data(tm_bg_file)
_, te_bg_wavenum, te_bg_single_beam, _ = load_data(te_bg_file)

old_meas = SinglePassMeas(ident=baseline_date)

old_meas.TE_single_beam=te_bg_single_beam
old_meas.TM_single_beam=tm_bg_single_beam
old_meas.TE_wavenum=te_bg_wavenum
old_meas.TM_wavenum=tm_bg_wavenum

axs[0].plot(old_meas.TE_wavenum, old_meas.TE_single_beam, label='TE ' + baseline_date,
            color='yellowgreen')
axs[0].plot(old_meas.TM_wavenum, old_meas.TM_single_beam, label='TM ' + baseline_date,
            color='darkgoldenrod')

blues = ['darkblue','mediumblue','blue','cornflowerblue']
reds = ['fuchsia','mediumvioletred','crimson','lightpink']


tm_file = os.path.join(base_dir, meas_date +'-bg-P0deg.CSV')
te_file = os.path.join(base_dir, meas_date+'-bg-P90deg.CSV')
_, tm_wavenum, tm_single_beam, _ = load_data(tm_file)
_, te_wavenum, te_single_beam, _ = load_data(te_file)
new_meas = SinglePassMeas(ident=meas_date)
new_meas.TM_wavenum = tm_wavenum
new_meas.TE_wavenum = te_wavenum
new_meas.TM_single_beam = tm_single_beam
new_meas.TE_single_beam = te_single_beam

#add the raw data plot
axs[0].plot(new_meas.TE_wavenum, new_meas.TE_single_beam, label='TE '+meas_date,
            color=blues[0])
axs[0].plot(new_meas.TM_wavenum, new_meas.TM_single_beam, label='TM '+meas_date,
            color=reds[0])

axs[1].plot(new_meas.TM_wavenum, new_meas.TM_single_beam/old_meas.TM_single_beam, label= 'TM',
            color=reds[0])
axs[1].plot(new_meas.TE_wavenum, new_meas.TE_single_beam/old_meas.TE_single_beam, label= 'TE',
            color=blues[0])

plt.figure(fig)
axs[0].legend()
axs[0].legend(prop={"size":14})
axs[1].legend()
axs[1].legend(prop={"size":14})

plt.tight_layout()
save_title = os.path.join(base_dir, baseline_date + ' vs ' + meas_date + 'raw scans and ratios' + '.svg')
plt.savefig(save_title)
#                    color='green')

plt.figure(fig)
plt.show()
