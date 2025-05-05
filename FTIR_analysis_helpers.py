import pandas as pd
# import matplotlib.pyplot as plt
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
        self.name = samp

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

class SinglePassMeas:
    def __init__(self, ident):
        self.TM_single_beam = None
        self.TE_single_beam = None
        self.TM_wavenum = None
        self.TE_wavenum = None

        self.TM_masked = None
        self.TE_masked = None
        self.TM_wavenum_masked = None
        self.TE_wavenum_masked = None

def build_MP(TEfile,TMfile,sample_name,nuextrema=None):
    _, tm_wavenum, tm_single_beam, _ = load_data(TMfile)
    _, te_wavenum, te_single_beam, _ = load_data(TEfile)
    samp_meas = MultipassMeas(samp=sample_name)

    # try using backgrounds with no sample in the path
    samp_meas.TM_wavenum = tm_wavenum
    samp_meas.TE_wavenum = te_wavenum
    samp_meas.TM_single_beam = tm_single_beam
    samp_meas.TE_single_beam = te_single_beam
    if nuextrema is not None:
        mask_samp = (samp_meas.TE_wavenum > nuextrema[0]) & (samp_meas.TE_wavenum < nuextrema[1])
        samp_meas.TE_masked = samp_meas.TE_single_beam[mask_samp]
        samp_meas.TM_masked = samp_meas.TM_single_beam[mask_samp]
        samp_meas.TM_wavenum_masked = samp_meas.TM_wavenum[mask_samp]
        samp_meas.TE_wavenum_masked = samp_meas.TE_wavenum[mask_samp]

    return samp_meas

def fitFnLorentz(nu, nuo, kappanuhalf, A, B):
    return A + ((B/np.pi) * kappanuhalf) / ((nu - nuo) ** 2 + kappanuhalf ** 2)

def fitFnNormal(nu, nuo, sigma, A, B):
    return A + (B/(sigma*np.sqrt(2*np.pi)))*np.exp((-nu+nuo)/(2*sigma**2))

def maskFit(nu_range,wavenum,alpha_ISB):
    mask_fit = (wavenum > nu_range[0]) & (wavenum < nu_range[1])
    alpha_ISB_select = alpha_ISB[mask_fit]
    wavenum_fit = wavenum[mask_fit]
    return alpha_ISB_select,wavenum_fit

def fitLorentzPlot(nu_range,kappanu_guess,wavenum,alpha_ISB,axs_fits,nu_fit_plot_range=None):

    alpha_ISB_select,wavenum_fit = maskFit(nu_range,wavenum,alpha_ISB)

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

    # fit_label = r"$\nu_0 = %0.2f {cm}^{-1}, \Delta \nu = %0.2f {cm}^{-1},A = %0.2f [units unknown], B= %0.2f $" % (fitLorentz[0], fitLorentz[1]*2,fitLorentz[2],fitLorentz[3])
    fit_label = r"$\nu_0 = %0.2f {cm}^{-1}$" % (fitLorentz[0]) + "\n" + r"$\Delta \nu = %0.2f {cm}^{-1}$" % (fitLorentz[1]*2)

    if nu_fit_plot_range is None:
        nu_fit_plot= wavenum
    else:
        trash,nu_fit_plot=maskFit(nu_fit_plot_range,wavenum,alpha_ISB)
    axs_fits.plot(nu_fit_plot, [fitFnLorentz(nu, *fitLorentz) for nu in nu_fit_plot],
                  linewidth=1,
                  label=fit_label)

    return fitLorentz

def fitNormalPlot(nu_range,mu_guess,sigma_guess,wavenum,alpha_ISB,axs_fits):

    alpha_ISB_select,wavenum_fit = maskFit(nu_range,wavenum,alpha_ISB)

    #calculate the fit guesses
    # nu0_guess = wavenum_fit[np.argmax(alpha_ISB_select)]
    A_guess = alpha_ISB_select[-1]
    B_guess = np.max(alpha_ISB_select) - np.min(alpha_ISB_select)
    fitGuess = (mu_guess, sigma_guess, A_guess,B_guess)

    nu0_bounds = [np.min(wavenum_fit),np.max(wavenum_fit)]

    sigma_bounds = [0,np.max(wavenum_fit) - np.min(wavenum_fit)]
    A_bounds = [-np.Infinity,np.Infinity]
    B_bounds = [0,np.Infinity]
    fitBounds = ([nu0_bounds[0], sigma_bounds[0], A_bounds[0], B_bounds[0]],
                 [nu0_bounds[1],sigma_bounds[1],A_bounds[1], B_bounds[1]])
    #
    fitNormal, trash = opt.curve_fit(fitFnNormal, wavenum_fit, alpha_ISB_select, p0=fitGuess,
                                      bounds=fitBounds)

    fit_label = r"$\mu = %0.2f {cm}^{-1}, \sigma = %0.2f {cm}^{-1},A = %0.2f [units unknown], B= %0.2f $" % (fitNormal[0], fitNormal[1]*2,fitNormal[2],fitNormal[3])
    axs_fits.plot(wavenum, [fitFnNormal(nu, *fitNormal) for nu in wavenum],
                  linewidth=1,
                  label=fit_label)

    return fitNormal

def calculate_Fresnels(theta_i_deg, n2, n1):
    # assumes theta_i in degrees
    theta_i = np.radians(theta_i_deg)
    sin_theta_t = n1 * np.sin(theta_i) / n2
    non_TIR_indices = np.where(sin_theta_t <= 1)
    # TIR_indices = np.where(sin_theta_t > 1)

    non_TIR_thetas = theta_i[non_TIR_indices]
    costhetat = np.sqrt(1 - sin_theta_t[non_TIR_indices] ** 2)

    numerator_p = n1 * costhetat - n2 * np.cos(non_TIR_thetas)
    denominator_p = n1 * costhetat + n2 * np.cos(non_TIR_thetas)

    Rp = np.full_like(theta_i,1.0)
    Rp_nonTIR = np.abs(numerator_p / denominator_p) ** 2
    Rp[non_TIR_indices] = Rp_nonTIR
    Tp = 1 - Rp

    numerator_s = n1 * np.cos(non_TIR_thetas) - n2 * costhetat
    denominator_s = n1 * np.cos(non_TIR_thetas) + n2 * costhetat
    Rs = np.full_like(theta_i,1.0)
    Rs_nonTIR = np.abs(numerator_s / denominator_s) ** 2
    Rs[non_TIR_indices] = Rs_nonTIR
    Ts = 1.0 - Rs


    # TIR_thetas = theta_i[TIR_indices]
    # print(Rp)

    return Rp, Tp, Rs, Ts, n1, n2
