import os
# import scipy
import matplotlib.pyplot as plt
import numpy as np
from FTIR_analysis_helpers import MultipassMeas
from FTIR_analysis_helpers import load_data
from FTIR_analysis_helpers import fitLorentzPlot
from FTIR_analysis_helpers import SinglePassMeas
# from ...FTIR_analysis_helpers import MultipassMeas
# from ...FTIR_analysis_helpers import fitLorentzPlot
from matplotlib.ticker import MaxNLocator
sample_name = 'P530'
# angle_of_incidence = 45

axislabelsfont=15
legendfont=15
ticksize=15
titlesize=20
linesize= 2

numax = 3100
numin = 785

#adjust with well thicknesses based on Lodo runsheet

base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250122_P530-A-SP'
combined_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250404_P530_combined'

angles = [230,240,250]
angle0 = 280

#raw data plots

fig_mc, axs_mc = plt.subplots(nrows=2,ncols=1,sharex=True,figsize=(8,10))

axs_mc[0].tick_params(axis='x',labelsize=ticksize)
axs_mc[0].tick_params(axis='y',labelsize=ticksize)
axs_mc[1].tick_params(axis='y',labelsize=ticksize)
axs_mc[1].tick_params(axis='x',labelsize=ticksize)

theta_variation_title =sample_name + ' Raw data'

fit_plot_title = sample_name+ " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax) + r" ${cm}^{-1}$"

fit_plot_title = sample_name+ " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax) + r" ${cm}^{-1} Samp vs. bg$"

tm_bg_file = os.path.join(base_dir, 'no_samp_P0deg' + '.CSV')
te_bg_file = os.path.join(base_dir, 'no_samp_P90deg' + '.CSV')

_, tm_bg_wavenum, tm_bg_single_beam, _ = load_data(tm_bg_file)
_, te_bg_wavenum, te_bg_single_beam, _ = load_data(te_bg_file)

bg_meas = SinglePassMeas(ident='background')

bg_meas.TE_single_beam=te_bg_single_beam
bg_meas.TM_single_beam=tm_bg_single_beam
bg_meas.TE_wavenum=te_bg_wavenum
bg_meas.TM_wavenum=tm_bg_wavenum

mask_samp = (bg_meas.TE_wavenum > numin) & (bg_meas.TE_wavenum < numax)

bg_meas.TE_masked = bg_meas.TE_single_beam[mask_samp]
bg_meas.TM_masked = bg_meas.TM_single_beam[mask_samp]
bg_meas.TM_wavenum_masked = bg_meas.TM_wavenum[mask_samp]
bg_meas.TE_wavenum_masked = bg_meas.TE_wavenum[mask_samp]


#absorption fit plot
fit_plot_title = "FTIR Absorption Spectra"

blues = ['darkblue','mediumblue','blue','cornflowerblue']
reds = ['darkmagenta','mediumorchid','deeppink','pink']

numins = [1042,1042,1042]
numaxs = [1272,1272,1272]
nuplot_range =[1000,1272]
kappa_nu_guesses = [30]
sim_transition_line_meV = 207.7
sim_transition_line = sim_transition_line_meV*8.065

for i in range(0,len(angles)):
    angle = angles[i]
    angle_filename = sample_name + '-SP-rot'+ str(angle)+'deg-'
    tm_file = os.path.join(base_dir, angle_filename+'P0deg' + '.CSV')
    te_file = os.path.join(base_dir, angle_filename+'P90deg' + '.CSV')
    _, tm_wavenum, tm_single_beam, _ = load_data(tm_file)
    _, te_wavenum, te_single_beam, _ = load_data(te_file)
    angle = angle0-angle
    angle_meas = SinglePassMeas(ident=str(angle))
    angle_meas.TM_wavenum = tm_wavenum
    angle_meas.TE_wavenum = te_wavenum
    angle_meas.TM_single_beam = tm_single_beam
    angle_meas.TE_single_beam = te_single_beam

    #add the raw data plot

    #do the masking
    angle_meas.TE_masked = angle_meas.TE_single_beam[mask_samp]
    angle_meas.TM_masked = angle_meas.TM_single_beam[mask_samp]
    angle_meas.TM_wavenum_masked = angle_meas.TM_wavenum[mask_samp]
    angle_meas.TE_wavenum_masked = angle_meas.TE_wavenum[mask_samp]

    #calculate the absorption coefficient times path length

    offset = np.log(bg_meas.TM_masked/bg_meas.TE_masked)
    # offset = bg_meas.TM_masked/bg_meas.TE_masked


    # alpha_ISB=-angle_meas.TM_masked/angle_meas.TE_masked + offset
    alpha_ISB= -np.log(angle_meas.TM_masked/angle_meas.TE_masked)+offset

    axs_mc[0].plot(angle_meas.TE_wavenum_masked, alpha_ISB, label= str(angle) + '$\degree$',
                            color=reds[i],linewidth=linesize)
    # closestidx = (np.abs(angle_meas.TE_wavenum_masked-sim_transition_line)).argmin()
    # if i ==0:
    #     axs_fits.axvline(sim_transition_line,label=r'$|d_{0,23}|=4.8 A$',color='lime')

    # do the fit
    nu_range = [numins[i],numaxs[i]]
    kappa_nu_guess = kappa_nu_guesses[0]
    # fitLorentzParams = fitLorentzPlot(nu_range, kappa_nu_guess,angle_meas.TE_wavenum_masked, alpha_ISB, axs_fits,nu_fit_plot_range=nuplot_range)

#get the 0 deg one in there
special_filename = 'P530-SP-top-chip'
special_file = os.path.join(base_dir, special_filename + '.CSV')
_, te_wavenum, te_single_beam, _ = load_data(special_file)
angle = 0
angle_meas = SinglePassMeas(ident=str(angle))
angle_meas.TE_wavenum = te_wavenum
angle_meas.TE_single_beam = te_single_beam

## add the P530F multipass ##
base_dir_P530F = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/P530F/20250307_P530F_MP'
sample_name_P530F = 'P530F'
Lsample = 6 #mm
sample_thick = 0.4215 #mm
Nbounces = Lsample/sample_thick
Lpath = Lsample*np.sqrt(2)
numax = 3100
numin = 1200

samp_meas = MultipassMeas(samp=sample_name_P530F)
bg_meas = MultipassMeas(samp='bg_ap_20_gain_4')

tm_file = os.path.join(base_dir_P530F, 'P530F_TM' + '.CSV')
te_file = os.path.join(base_dir_P530F, 'P530F_TE' + '.CSV')

te_bg_file = os.path.join(base_dir_P530F,'bg_TE_ap_20_gain_4'+'.CSV')
tm_bg_file = os.path.join(base_dir_P530F,'bg_TM_ap_20_gain_4'+'.CSV')

_, tm_wavenum, tm_single_beam, _ = load_data(tm_file)
_, te_wavenum, te_single_beam, _ = load_data(te_file)

_, tm_bg_wavenum, tm_bg_single_beam, _ = load_data(tm_bg_file)
_, te_bg_wavenum, te_bg_single_beam, _ = load_data(te_bg_file)

samp_meas.TM_wavenum = tm_wavenum
samp_meas.TE_wavenum = te_wavenum
samp_meas.TM_single_beam = tm_single_beam
samp_meas.TE_single_beam = te_single_beam

bg_meas.TM_wavenum = tm_bg_wavenum
bg_meas.TE_wavenum = te_bg_wavenum
bg_meas.TM_single_beam = tm_bg_single_beam
bg_meas.TE_single_beam = te_bg_single_beam

mask_samp = (samp_meas.TE_wavenum > numin) & (samp_meas.TE_wavenum < numax)
mask_bg = (bg_meas.TE_wavenum > numin) & (bg_meas.TE_wavenum < numax)

samp_meas.TE_masked = samp_meas.TE_single_beam[mask_samp]
samp_meas.TM_masked = samp_meas.TM_single_beam[mask_samp]
samp_meas.TM_wavenum_masked = samp_meas.TM_wavenum[mask_samp]
samp_meas.TE_wavenum_masked = samp_meas.TE_wavenum[mask_samp]

bg_meas.TE_masked = bg_meas.TE_single_beam[mask_bg]
bg_meas.TM_masked = bg_meas.TM_single_beam[mask_bg]
bg_meas.TM_wavenum_masked = bg_meas.TM_wavenum[mask_bg]
bg_meas.TE_wavenum_masked = bg_meas.TE_wavenum[mask_bg]

offset = np.log(bg_meas.TM_masked / bg_meas.TE_masked)

alpha_ISB = -np.log(samp_meas.TM_masked / samp_meas.TE_masked) + offset

path_len_label = r"$N_{bounces} = %d $" % (Nbounces) + r"$, L_{path} = %0.2f mm$" % (Lpath)

axs_mc[1].plot(samp_meas.TE_wavenum_masked, alpha_ISB, label=path_len_label,
              color='b',linewidth=linesize)

## add the straight on P530-C single pass ##

sample_name = 'P530-C'
numax = 3100
numin = 800
ap = 20
gain = 2
Fout = 0.78 #fresnel coefficient in
Fin = Fout
nu_ISB = 1120

angle=0

settings_suffix = 'ap-'+str(ap)+'-gain-'+str(gain)
base_dir_P530C = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250204_P530-C-SP'


tm_bg_file = os.path.join(base_dir_P530C, 'bg-P0deg-' + settings_suffix + '.CSV')
te_bg_file = os.path.join(base_dir_P530C, 'bg-P90deg-' + settings_suffix + '.CSV')

_, tm_bg_wavenum, tm_bg_single_beam, _ = load_data(tm_bg_file)
_, te_bg_wavenum, te_bg_single_beam, _ = load_data(te_bg_file)

bg_meas = SinglePassMeas(ident='background')

bg_meas.TE_single_beam=te_bg_single_beam
bg_meas.TM_single_beam=tm_bg_single_beam
bg_meas.TE_wavenum=te_bg_wavenum
bg_meas.TM_wavenum=tm_bg_wavenum

mask_samp = (bg_meas.TE_wavenum > numin) & (bg_meas.TE_wavenum < numax)

bg_meas.TE_masked = bg_meas.TE_single_beam[mask_samp]
bg_meas.TM_masked = bg_meas.TM_single_beam[mask_samp]
bg_meas.TM_wavenum_masked = bg_meas.TM_wavenum[mask_samp]
bg_meas.TE_wavenum_masked = bg_meas.TE_wavenum[mask_samp]

angle_filename = sample_name + '-SP-rot-' + str(angle) + 'deg'

tm_file = os.path.join(base_dir_P530C, angle_filename + '-P0deg-' + settings_suffix + '.CSV')
te_file = os.path.join(base_dir_P530C, angle_filename + '-P90deg-' + settings_suffix + '.CSV')
_, tm_wavenum, tm_single_beam, _ = load_data(tm_file)
_, te_wavenum, te_single_beam, _ = load_data(te_file)
angle_meas = SinglePassMeas(ident=str(angle))
angle_meas.TM_wavenum = tm_wavenum
angle_meas.TE_wavenum = te_wavenum
angle_meas.TM_single_beam = tm_single_beam
angle_meas.TE_single_beam = te_single_beam

angle_meas.TE_masked = angle_meas.TE_single_beam[mask_samp]
angle_meas.TM_masked = angle_meas.TM_single_beam[mask_samp]
angle_meas.TM_wavenum_masked = angle_meas.TM_wavenum[mask_samp]
angle_meas.TE_wavenum_masked = angle_meas.TE_wavenum[mask_samp]

offset = np.log(bg_meas.TM_masked / bg_meas.TE_masked)

alpha_ISB = -np.log(angle_meas.TM_masked / angle_meas.TE_masked) + offset
axs_mc[0].plot(angle_meas.TE_wavenum_masked, alpha_ISB, label= str(angle) + '$\degree$',color=reds[3],linewidth=linesize)

## do overall figure formatting ##

plt.figure(fig_mc)
fig_mc.supxlabel(r"Wavenumber [${cm}^{-1}$]",fontsize=axislabelsfont)
# fig_mc.supylabel(r"$\alpha_{ISB} \times L_{path}=-\ln (\frac{I_{out,TM}}{I_{out,TE}}) + \ln(\frac{I_{bg,TM}}{I_{bg,TE}})$",fontsize=axislabelsfont)
fig_mc.suptitle('P530 Absorption: Single Pass and Multipass',fontsize=titlesize)
axs_mc[0].legend(loc='upper right')
axs_mc[1].legend(loc='upper right')
# axs_mc[2].legend(loc='upper right')
axs_mc[0].legend(prop={"size":legendfont})
axs_mc[1].legend(prop={"size":legendfont})
# axs_mc[2].legend(prop={"size":legendfont})
axs_mc[0].grid()
axs_mc[1].grid()
# axs_mc[2].grid()

axs_mc[0].xaxis.set_major_locator(MaxNLocator(integer=True))

plt.tight_layout()
save_title_mc = os.path.join(combined_dir, sample_name + 'mc_combined' + '.svg')
plt.savefig(save_title_mc)
plt.show()