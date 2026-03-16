import pandas as pd
# import matplotlib.pyplot as plt
import os
import numpy as np
import scipy.optimize as opt
import re
from datetime import datetime
from dataclasses import dataclass, field

def load_data(filename, return_date=False):
    data = pd.read_csv(filename, header=None)
    wavenumber = np.array(data[0], dtype=float)
    single_beam = np.array(data[1], dtype=float)
    wavelength = 1e4 / wavenumber  # cm^-1 -> µm

    label = os.path.splitext(os.path.basename(filename))[0]

    parsed_date = None

    if return_date:
        # Find all 8-digit sequences in the full path
        matches = re.findall(r'\d{8}', filename)

        if matches:
            date_str = matches[-1]  # right-most wins
            try:
                parsed_date = datetime.strptime(date_str, "%Y%m%d").date()
            except ValueError:
                parsed_date = None
        else:
            parsed_date = None

        return wavelength, wavenumber, single_beam, label, parsed_date
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
    def __init__(self, samp,thetai):
        self.name = samp
        self.thetai = thetai

        self.TM_single_beam = None
        self.TE_single_beam = None
        self.TM_wavenum = None
        self.TE_wavenum = None

        self.TM_masked = None
        self.TE_masked = None
        self.TM_wavenum_masked = None
        self.TE_wavenum_masked = None

def calculate_Fresnels(theta_i_deg, n1, n2):
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

    # theta_t = np.full_like(theta_i,90.)
    theta_t_nonTIR = np.rad2deg(np.arcsin(sin_theta_t[non_TIR_indices]))
    # theta_t[non_TIR_indices] = theta_t_nonTIR


    # TIR_thetas = theta_i[TIR_indices]
    # print(Rp)

    return Rp, Tp, Rs, Ts, n1, n2,theta_t_nonTIR

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

def build_SP(TEfile,TMfile,sample_name,thetai,fresnel=False,n1=None,n2=None,n3=None,nuextrema=None):
    _, tm_wavenum, tm_single_beam, _ = load_data(TMfile)
    _, te_wavenum, te_single_beam, _ = load_data(TEfile)
    samp_meas = SinglePassMeas(samp=sample_name,thetai=thetai)

    # try using backgrounds with no sample in the path
    samp_meas.TM_wavenum = tm_wavenum
    samp_meas.TE_wavenum = te_wavenum
    samp_meas.TM_single_beam_raw = tm_single_beam
    samp_meas.TE_single_beam_raw = te_single_beam

    #account for fresnel
    if fresnel:
        #calc fresnel coefficients for beam from air to sampe
        _, Tp12, _, Ts12,_ , _,theta_t_nonTIR12 = calculate_Fresnels(thetai,n1,n2)
        print("TP12: " + str(Tp12) + " Ts12: " + str(Ts12))

        #calcualte fresnel for samp to air
        _, Tp23, _, Ts23, _, _, theta_t_nonTIR23 = calculate_Fresnels(theta_t_nonTIR12, n2, n3)
        print("TP23: " + str(Tp23) + " Ts23: " + str(Ts23))
        #apply to relevant angles

        samp_meas.TM_single_beam = samp_meas.TM_single_beam_raw/(Tp12*Tp23)
        samp_meas.TE_single_beam = samp_meas.TE_single_beam_raw/(Ts12*Ts23)
    else:
        samp_meas.TM_single_beam = samp_meas.TM_single_beam_raw
        samp_meas.TE_single_beam = samp_meas.TE_single_beam_raw

    if nuextrema is not None:
        mask_samp = (samp_meas.TE_wavenum > nuextrema[0]) & (samp_meas.TE_wavenum < nuextrema[1])
        samp_meas.TE_masked = samp_meas.TE_single_beam[mask_samp]
        samp_meas.TM_masked = samp_meas.TM_single_beam[mask_samp]
        samp_meas.TM_wavenum_masked = samp_meas.TM_wavenum[mask_samp]
        samp_meas.TE_wavenum_masked = samp_meas.TE_wavenum[mask_samp]

    return samp_meas

@dataclass
class FitResult:
    nu0s:    np.ndarray  # shape (N,) — fitted center wavenumbers
    kappas:  np.ndarray  # shape (N,) — fitted HWHM for each peak
    A:       float       # fitted baseline offset
    Bs:      np.ndarray  # shape (N,) — fitted amplitudes
    cov:     np.ndarray  # covariance matrix from curve_fit
    n_peaks: int
    nu_range: list       # [nu_min, nu_max] window used for fitting


def fitFnLorentz(nu, nuo, kappanuhalf, A, B):
    return A + ((B/np.pi) * kappanuhalf) / ((nu - nuo) ** 2 + kappanuhalf ** 2)


def fitFnNLorentz(nu, *params):
    # params: [nu0_1, kappa_1, B_1,  nu0_2, kappa_2, B_2,  ...,  A]
    A = params[-1]
    n_peaks = (len(params) - 1) // 3
    result = A
    for i in range(n_peaks):
        nu0, kappa, B = params[3*i], params[3*i + 1], params[3*i + 2]
        result = result + (B/np.pi * kappa) / ((nu - nu0)**2 + kappa**2)
    return result

def fitFnNormal(nu, nuo, sigma, A, B):
    return A + (B/(sigma*np.sqrt(2*np.pi)))*np.exp((-nu+nuo)/(2*sigma**2))

def maskFit(nu_range, wavenum, alpha_ISB):
    mask_fit = (wavenum > nu_range[0]) & (wavenum < nu_range[1])
    alpha_ISB_select = alpha_ISB[mask_fit]
    wavenum_fit = wavenum[mask_fit]
    return alpha_ISB_select, wavenum_fit


def fit_lorentz(nu_range, nu0_guess, kappa_guess, wavenum, alpha_ISB,
                A_guess=None, B_guess=None):
    alpha_ISB_select, wavenum_fit = maskFit(nu_range, wavenum, alpha_ISB)

    if A_guess is None:
        A_guess = alpha_ISB_select[-1]
    if B_guess is None:
        B_guess = np.max(alpha_ISB_select) - np.min(alpha_ISB_select)

    p0 = [nu0_guess, kappa_guess, A_guess, B_guess]
    nu_width = nu_range[1] - nu_range[0]
    bounds = ([nu_range[0], 0,       -np.inf, 0      ],
              [nu_range[1], nu_width, np.inf,  np.inf ])

    params, cov = opt.curve_fit(fitFnLorentz, wavenum_fit, alpha_ISB_select,
                                p0=p0, bounds=bounds)
    return FitResult(nu0s=np.array([params[0]]), kappas=np.array([params[1]]),
                     A=params[2], Bs=np.array([params[3]]),
                     cov=cov, n_peaks=1, nu_range=nu_range)


def fit_nlorentz(nu_range, nu0_guesses, kappa_guesses, wavenum, alpha_ISB,
                 A_guess=None, B_guesses=None):
    alpha_ISB_select, wavenum_fit = maskFit(nu_range, wavenum, alpha_ISB)
    n = len(nu0_guesses)

    if A_guess is None:
        A_guess = alpha_ISB_select[-1]
    if B_guesses is None:
        peak_height = np.max(alpha_ISB_select) - np.min(alpha_ISB_select)
        B_guesses = [peak_height / n] * n

    # pack as [nu0_1, kappa_1, B_1,  nu0_2, kappa_2, B_2,  ...,  A]
    p0 = []
    for i in range(n):
        p0 += [nu0_guesses[i], kappa_guesses[i], B_guesses[i]]
    p0.append(A_guess)

    nu_width = nu_range[1] - nu_range[0]
    lower, upper = [], []
    for _ in range(n):
        lower += [nu_range[0], 0,        0      ]
        upper += [nu_range[1], nu_width, np.inf ]
    lower.append(-np.inf)
    upper.append(np.inf)

    params, cov = opt.curve_fit(fitFnNLorentz, wavenum_fit, alpha_ISB_select,
                                p0=p0, bounds=(lower, upper))

    nu0s   = np.array([params[3*i]     for i in range(n)])
    kappas = np.array([params[3*i + 1] for i in range(n)])
    Bs     = np.array([params[3*i + 2] for i in range(n)])
    return FitResult(nu0s=nu0s, kappas=kappas, A=params[-1], Bs=Bs,
                     cov=cov, n_peaks=n, nu_range=nu_range)


def plot_lorentz_fit(fit_result, wavenum, axs, nu_fit_plot_range=None):
    if nu_fit_plot_range is None:
        nu_plot = wavenum
    else:
        mask = (wavenum > nu_fit_plot_range[0]) & (wavenum < nu_fit_plot_range[1])
        nu_plot = wavenum[mask]

    # build full combined params: [nu0_1, kappa_1, B_1, ..., A]
    params = []
    for i in range(fit_result.n_peaks):
        params += [fit_result.nu0s[i], fit_result.kappas[i], fit_result.Bs[i]]
    params.append(fit_result.A)

    # for N > 1: plot each individual component (baseline + single peak) as dashed
    if fit_result.n_peaks > 1:
        for i in range(fit_result.n_peaks):
            nu0, kappa = fit_result.nu0s[i], fit_result.kappas[i]
            component_params = [nu0, kappa, fit_result.A, fit_result.Bs[i]]
            peak_label = (
                r"$\nu_0 = %.2f\ \mathrm{cm}^{-1}$" % nu0 + "\n"
                + r"$\Delta\nu = %.2f\ \mathrm{cm}^{-1}$" % (kappa * 2) + "\n"
                + r"$\Delta\nu/\nu = %.2f\%%$" % (kappa / nu0 * 2 * 100)
            )
            axs.plot(nu_plot, [fitFnLorentz(nu, *component_params) for nu in nu_plot],
                     linewidth=1, linestyle='--', label=peak_label)
        combined_label = 'combined fit'
    else:
        nu0, kappa = fit_result.nu0s[0], fit_result.kappas[0]
        combined_label = (
            r"$\nu_0 = %.2f\ \mathrm{cm}^{-1}$" % nu0 + "\n"
            + r"$\Delta\nu = %.2f\ \mathrm{cm}^{-1}$" % (kappa * 2) + "\n"
            + r"$\Delta\nu/\nu = %.2f\%%$" % (kappa / nu0 * 2 * 100)
        )

    axs.plot(nu_plot, [fitFnNLorentz(nu, *params) for nu in nu_plot],
             linewidth=1, label=combined_label)


def fitLorentzPlot(nu_range, kappanu_guess, wavenum, alpha_ISB, axs_fits,
                   nu_fit_plot_range=None):
    """Backwards-compatible wrapper around fit_lorentz + plot_lorentz_fit."""
    alpha_ISB_select, wavenum_fit = maskFit(nu_range, wavenum, alpha_ISB)
    nu0_guess = wavenum_fit[np.argmax(alpha_ISB_select)]
    result = fit_lorentz(nu_range, nu0_guess, kappanu_guess, wavenum, alpha_ISB)
    plot_lorentz_fit(result, wavenum, axs_fits, nu_fit_plot_range)
    return np.array([result.nu0s[0], result.kappas[0], result.A, result.Bs[0]])

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

