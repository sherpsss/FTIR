import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import scipy.optimize as opt

def load_data(filename):
    data = pd.read_csv(filename, header=None)
    wavenumber = np.array(data[0])# First column: wavenumber
    single_beam = np.array(data[1])  # Second column: transmission
    wavelength = 1e4 / wavenumber  # Convert wavenumber (cm^-1) to wavelength (µm)

    # Extract a label from the filename (remove path and extension)
    label = os.path.splitext(os.path.basename(filename))[0]

    return wavelength, wavenumber, single_beam, label

class MultipassMeas:
    def __init__(self,samp):
        self.TM_single_beam = None
        self.TE_single_beam = None
        self.TM_wavenum = None
        self.TE_wavenum = None

        self.TM_reshaped = None
        self.TE_reshaped = None
        self.TE_wavenum_reshaped = None
        self.TM_wavenum_reshaped = None

        self.TM_masked = None
        self.TE_masked = None
        self.TM_wavenum_masked = None
        self.TE_wavenum_masked = None


def fitFnLorentz(nu, nuo, kappanuhalf, A, B):

    return A + ((B/np.pi) * kappanuhalf) / ((nu - nuo) ** 2 + kappanuhalf ** 2)

def fitLorentzPlot(nu_range,kappanu_guess,wavenum,alpha_ISB,axs_fits):

    mask_fit = (wavenum > nu_range[0]) & (wavenum < nu_range[1])
    alpha_ISB_select = alpha_ISB[mask_fit]
    wavenum_fit = wavenum[mask_fit]

    #calculate the fit guesses
    nu0_guess = wavenum_fit[np.argmax(alpha_ISB_select)]
    A_guess = alpha_ISB_select[-1]
    B_guess = np.max(alpha_ISB_select) - np.min(alpha_ISB_select)
    fitGuess = (nu0_guess, kappanu_guess, A_guess,B_guess)

    nu0_bounds = [np.min(wavenum_fit),np.max(wavenum_fit)]
    # kappanu_bounds = [0,np.max(wavenum_fit) - np.min(wavenum_fit)]
    # A_bounds = [-2 * np.abs(np.min(alpha_ISB)),2 * np.abs(np.min(alpha_ISB))]
    # B_bounds = [np.max(np.abs(alpha_ISB_select)) / 4,np.max(np.abs(alpha_ISB_select)) * 2]
    # fitBounds = ([nu0_bounds[0], kappanu_bounds[0], A_bounds[0], B_bounds[0]],
    #              [nu0_bounds[1],kappanu_bounds[1],A_bounds[1], B_bounds[1]])

    kappanu_bounds = [0,np.max(wavenum_fit) - np.min(wavenum_fit)]
    A_bounds = [-np.Infinity,np.Infinity]
    B_bounds = [0,np.Infinity]
    fitBounds = ([nu0_bounds[0], kappanu_bounds[0], A_bounds[0], B_bounds[0]],
                 [nu0_bounds[1],kappanu_bounds[1],A_bounds[1], B_bounds[1]])
    #
    fitLorentz, trash = opt.curve_fit(fitFnLorentz, wavenum_fit, alpha_ISB_select, p0=fitGuess,
                                      bounds=fitBounds)

    fit_label = r"$\nu_0 = %0.2f {cm}^{-1}, \Delta \nu = %0.2f {cm}^{-1},A = %0.2f [units unknown], B= %0.2f $" % (fitLorentz[0], fitLorentz[1]*2,fitLorentz[2],fitLorentz[3])
    axs_fits.plot(wavenum, [fitFnLorentz(nu, *fitLorentz) for nu in wavenum],
                  linewidth=1,
                  label=fit_label)

    return fitLorentz


sample_name = 'P12-3-24-1'
Lsample = 10 #mm
# angle_of_incidence = 45
sample_thick = 0.5 #mm
Nbounces = Lsample/sample_thick
well_period = 320e-5 #mm
Nperiods = 29
epi_path = np.sqrt(2)*Nperiods*well_period
Lpath = epi_path*Nbounces

numax = 3100
numin = 750

#adjust with well thicknesses based on Lodo runsheet

samp_meas = MultipassMeas(samp=sample_name)

base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/P12-3-24-1'

tm_file = os.path.join(base_dir, 'P0deg' + '.CSV')
te_file = os.path.join(base_dir, 'P90deg' + '.CSV')

_, tm_wavenum, tm_single_beam, _ = load_data(tm_file)
_, te_wavenum, te_single_beam, _ = load_data(te_file)

#try using backgrounds with no sample in the path
samp_meas.TM_wavenum = tm_wavenum
samp_meas.TE_wavenum = te_wavenum
samp_meas.TM_single_beam = tm_single_beam
samp_meas.TE_single_beam = te_single_beam

background_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20241211_theta_i_sweep_singlepass_342A/backgrounds'
bg_name = 'no_samp'
bg_meas = MultipassMeas(samp=bg_name)

tm_bg_file = os.path.join(background_dir, 'background_P0deg' + '.CSV')
te_bg_file = os.path.join(background_dir, 'background_P90deg' + '.CSV')

_, tm_bg_wavenum, tm_bg_single_beam, _ = load_data(tm_bg_file)
_, te_bg_wavenum, te_bg_single_beam, _ = load_data(te_bg_file)

bg_meas.TE_single_beam=te_bg_single_beam
bg_meas.TM_single_beam=tm_bg_single_beam
bg_meas.TE_wavenum=te_bg_wavenum
bg_meas.TM_wavenum=tm_bg_wavenum

fig, axs = plt.subplots(1, 2, figsize=(10, 10))

axs[0].plot(samp_meas.TE_wavenum, samp_meas.TE_single_beam, label= 'TE',
                        color='blue')
axs[0].plot(samp_meas.TM_wavenum , samp_meas.TM_single_beam , label= 'TM',
                        color='red')

axs[0].plot(bg_meas.TE_wavenum, bg_meas.TE_single_beam, label= 'TE bg ' + bg_name,
                        color='c')
axs[0].plot(bg_meas.TM_wavenum , bg_meas.TM_single_beam , label= 'TM bg' + bg_name,
                        color='m')

samp_meas.TE_reshaped = samp_meas.TE_single_beam.reshape(-1,2).mean(axis=1)
samp_meas.TM_reshaped = samp_meas.TM_single_beam.reshape(-1,2).mean(axis=1)
samp_meas.TE_wavenum_reshaped = samp_meas.TE_wavenum.reshape(-1,2).mean(axis=1)
samp_meas.TM_wavenum_reshaped = samp_meas.TM_wavenum.reshape(-1,2).mean(axis=1)

axs[0].plot(bg_meas.TE_wavenum, samp_meas.TE_reshaped, label= 'TE  reshaped',
                        color='dodgerblue')
axs[0].plot(bg_meas.TM_wavenum , samp_meas.TM_reshaped , label= 'TM  reshaped',
                        color='sienna')

#now do the masking

mask_samp = (samp_meas.TE_wavenum_reshaped > numin) & (samp_meas.TE_wavenum_reshaped < numax)
#
# # Apply the mask to the array
#
mask_bg = (bg_meas.TE_wavenum > numin) & (bg_meas.TE_wavenum < numax)



#average the signal to compare


axs[0].set_xlabel('Wavenumber (cm^-1)')
axs[0].set_ylabel("Single Beam")
theta_variation_title =sample_name
axs[0].set_title(theta_variation_title)

fit_plot_title = sample_name+ " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax) + r" ${cm}^{-1}$"
axs[1].set_title(fit_plot_title)
axs[1].set_xlabel('Wavenumber (cm^-1)')
axs[1].set_ylabel('Transmission Ratio')
theta_variation_title_polarization_ratios = theta_variation_title + ' polarization ratios'

samp_meas.TE_masked = samp_meas.TE_reshaped[mask_samp]
samp_meas.TM_masked = samp_meas.TM_reshaped[mask_samp]
samp_meas.TM_wavenum_masked = samp_meas.TM_wavenum_reshaped[mask_samp]
samp_meas.TE_wavenum_masked = samp_meas.TE_wavenum_reshaped[mask_samp]

axs[1].plot(samp_meas.TE_wavenum_masked, samp_meas.TM_masked/samp_meas.TE_masked, label= 'TM/TE with sample masked reshaped',
                        color='green')
axs[1].plot(bg_meas.TE_wavenum[mask_bg], bg_meas.TM_single_beam[mask_bg]/bg_meas.TE_single_beam[mask_bg], label= 'TM/TE no sample masked reshaped',
                        color='y')

axs[0].legend()
axs[1].legend()
plt.tight_layout()
plt.show()

#fit figure

fig_fits, axs_fits = plt.subplots(figsize=(10, 10))
axs_fits.set_xlabel('Wavenumber (cm^-1)')

offset = np.log(bg_meas.TM_single_beam[mask_bg]/bg_meas.TE_single_beam[mask_bg])

alpha_ISB= -np.log(samp_meas.TM_masked/samp_meas.TE_masked)+offset

axs_fits.plot(samp_meas.TE_wavenum_masked, alpha_ISB, label= r"$-\ln (\frac{I_{out,TM}}{I_{out,TE}}) + \ln(\frac{I_{bg,TM}}{I_{bg,TE}})$",
                        color='green')

# #do fits
# numins = [840,1012,1700,1881]
# numaxs = [1039,1700,1780,1910]
numins = [860,1012,1700,1881]
numaxs = [1039,1700,1780,1910]
kappa_nu_guesses = [35,50,10,5]
# numins = [840,1012]
# numaxs = [1039,1700]
for i in range(0,len(numins)):
    nu_range = [numins[i],numaxs[i]]
    kappa_nu_guess = kappa_nu_guesses[i]
    fitLorentzPlot(nu_range, kappa_nu_guess,samp_meas.TE_wavenum_masked, alpha_ISB, axs_fits)


axs_fits.set_ylabel(r"$\alpha_{ISB} \times L_{path}$ ")
fit_plot_title = sample_name+ " SNR mask " + str(numin) + r"$ < \nu < $" + str(numax)
axs_fits.set_title(fit_plot_title)
axs_fits.legend()
plt.figure(fig_fits)
plt.tight_layout()
save_title = os.path.join(base_dir, sample_name + 'Lorentzian fits' + '.svg')
plt.savefig(save_title)
plt.show()
