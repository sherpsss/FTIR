import pandas as pd
import os
import numpy as np
import scipy.optimize as opt
import re
from datetime import datetime
from dataclasses import dataclass


# ── data loading ──────────────────────────────────────────────────────────────

def load_data(filename, return_date=False):
    data = pd.read_csv(filename, header=None)
    wavenumber   = np.array(data[0], dtype=float)
    single_beam  = np.array(data[1], dtype=float)
    wavelength   = 1e4 / wavenumber  # cm^-1 -> µm
    label        = os.path.splitext(os.path.basename(filename))[0]

    if return_date:
        matches = re.findall(r'\d{8}', filename)
        parsed_date = None
        if matches:
            try:
                parsed_date = datetime.strptime(matches[-1], "%Y%m%d").date()
            except ValueError:
                pass
        return wavelength, wavenumber, single_beam, label, parsed_date

    return wavelength, wavenumber, single_beam, label


# ── fit result container ──────────────────────────────────────────────────────

@dataclass
class FitResult:
    nu0s:     np.ndarray   # shape (N,) — fitted center wavenumbers
    kappas:   np.ndarray   # shape (N,) — fitted HWHM for each peak
    A:        float        # fitted baseline offset
    Bs:       np.ndarray   # shape (N,) — fitted amplitudes
    cov:      np.ndarray   # covariance matrix from curve_fit
    n_peaks:  int
    nu_range: list         # [nu_min, nu_max] window used for fitting


# ── line shape functions ──────────────────────────────────────────────────────

def fitFnLorentz(nu, nuo, kappanuhalf, A, B):
    return A + ((B / np.pi) * kappanuhalf) / ((nu - nuo) ** 2 + kappanuhalf ** 2)


def fitFnNLorentz(nu, *params):
    # params layout: [nu0_1, kappa_1, B_1,  nu0_2, kappa_2, B_2,  ...,  A]
    A       = params[-1]
    n_peaks = (len(params) - 1) // 3
    result  = A
    for i in range(n_peaks):
        nu0, kappa, B = params[3*i], params[3*i + 1], params[3*i + 2]
        result += (B / np.pi * kappa) / ((nu - nu0) ** 2 + kappa ** 2)
    return result


def fitFnNormal(nu, nuo, sigma, A, B):
    return A + (B / (sigma * np.sqrt(2 * np.pi))) * np.exp((-nu + nuo) / (2 * sigma ** 2))


# ── fitting helpers ───────────────────────────────────────────────────────────

def maskFit(nu_range, wavenum, alpha_ISB):
    mask = (wavenum > nu_range[0]) & (wavenum < nu_range[1])
    return alpha_ISB[mask], wavenum[mask]


def fit_lorentz(nu_range, nu0_guess, kappa_guess, wavenum, alpha_ISB,
                A_guess=None, B_guess=None):
    alpha_select, wavenum_fit = maskFit(nu_range, wavenum, alpha_ISB)

    if A_guess is None:
        A_guess = alpha_select[-1]
    if B_guess is None:
        B_guess = np.max(alpha_select) - np.min(alpha_select)

    p0       = [nu0_guess, kappa_guess, A_guess, B_guess]
    nu_width = nu_range[1] - nu_range[0]
    bounds   = ([nu_range[0], 0,        -np.inf, 0      ],
                [nu_range[1], nu_width,  np.inf, np.inf ])

    params, cov = opt.curve_fit(fitFnLorentz, wavenum_fit, alpha_select,
                                p0=p0, bounds=bounds)
    return FitResult(nu0s=np.array([params[0]]), kappas=np.array([params[1]]),
                     A=params[2], Bs=np.array([params[3]]),
                     cov=cov, n_peaks=1, nu_range=nu_range)


def fit_nlorentz(nu_range, nu0_guesses, kappa_guesses, wavenum, alpha_ISB,
                 A_guess=None, B_guesses=None):
    alpha_select, wavenum_fit = maskFit(nu_range, wavenum, alpha_ISB)
    n = len(nu0_guesses)

    if A_guess is None:
        A_guess = alpha_select[-1]
    if B_guesses is None:
        B_guesses = [(np.max(alpha_select) - np.min(alpha_select)) / n] * n

    p0 = []
    for i in range(n):
        p0 += [nu0_guesses[i], kappa_guesses[i], B_guesses[i]]
    p0.append(A_guess)

    nu_width     = nu_range[1] - nu_range[0]
    lower, upper = [], []
    for _ in range(n):
        lower += [nu_range[0], 0,        0      ]
        upper += [nu_range[1], nu_width, np.inf ]
    lower.append(-np.inf)
    upper.append(np.inf)

    params, cov = opt.curve_fit(fitFnNLorentz, wavenum_fit, alpha_select,
                                p0=p0, bounds=(lower, upper))

    nu0s   = np.array([params[3*i]     for i in range(n)])
    kappas = np.array([params[3*i + 1] for i in range(n)])
    Bs     = np.array([params[3*i + 2] for i in range(n)])
    return FitResult(nu0s=nu0s, kappas=kappas, A=params[-1], Bs=Bs,
                     cov=cov, n_peaks=n, nu_range=nu_range)


# ── plot helper ───────────────────────────────────────────────────────────────

def plot_lorentz_fit(fit_result, wavenum, axs, nu_fit_plot_range=None, lw=1.0,
                     component_colors=None, combined_color=None):
    if nu_fit_plot_range is None:
        nu_plot = wavenum
    else:
        mask    = (wavenum > nu_fit_plot_range[0]) & (wavenum < nu_fit_plot_range[1])
        nu_plot = wavenum[mask]

    # full combined params: [nu0_1, kappa_1, B_1, ..., A]
    params = []
    for i in range(fit_result.n_peaks):
        params += [fit_result.nu0s[i], fit_result.kappas[i], fit_result.Bs[i]]
    params.append(fit_result.A)

    if fit_result.n_peaks > 1:
        for i in range(fit_result.n_peaks):
            nu0, kappa       = fit_result.nu0s[i], fit_result.kappas[i]
            component_params = [nu0, kappa, fit_result.A, fit_result.Bs[i]]
            peak_label       = (
                r"$\nu_0 = %.1f\ \mathrm{cm}^{-1}$" % nu0 + "\n"
                + r"$\Delta\nu = %.1f\ \mathrm{cm}^{-1}$" % (kappa * 2) + "\n"
                + r"$\Delta\nu/\nu = %.1f\%%$" % (kappa / nu0 * 2 * 100)
            )
            color_kw = {'color': component_colors[i]} if component_colors is not None else {}
            axs.plot(nu_plot, [fitFnLorentz(nu, *component_params) for nu in nu_plot],
                     linewidth=lw, linestyle='--', label=peak_label, **color_kw)
        combined_label = 'combined fit'
    else:
        nu0, kappa     = fit_result.nu0s[0], fit_result.kappas[0]
        combined_label = (
            r"$\nu_0 = %.1f\ \mathrm{cm}^{-1}$" % nu0 + "\n"
            + r"$\Delta\nu = %.1f\ \mathrm{cm}^{-1}$" % (kappa * 2) + "\n"
            + r"$\Delta\nu/\nu = %.1f\%%$" % (kappa / nu0 * 2 * 100)
        )

    combined_color_kw = {'color': combined_color} if combined_color is not None else {}
    axs.plot(nu_plot, [fitFnNLorentz(nu, *params) for nu in nu_plot],
             linewidth=lw, label=combined_label, linestyle='--', **combined_color_kw)


# ── Fresnel coefficients ──────────────────────────────────────────────────────

def calculate_Fresnels(theta_i_deg, n1, n2):
    theta_i      = np.radians(theta_i_deg)
    sin_theta_t  = n1 * np.sin(theta_i) / n2
    non_TIR      = np.where(sin_theta_t <= 1)
    non_TIR_thetas = theta_i[non_TIR]
    costhetat    = np.sqrt(1 - sin_theta_t[non_TIR] ** 2)

    num_p = n1 * costhetat - n2 * np.cos(non_TIR_thetas)
    den_p = n1 * costhetat + n2 * np.cos(non_TIR_thetas)
    Rp    = np.full_like(theta_i, 1.0)
    Rp[non_TIR] = np.abs(num_p / den_p) ** 2
    Tp    = 1 - Rp

    num_s = n1 * np.cos(non_TIR_thetas) - n2 * costhetat
    den_s = n1 * np.cos(non_TIR_thetas) + n2 * costhetat
    Rs    = np.full_like(theta_i, 1.0)
    Rs[non_TIR] = np.abs(num_s / den_s) ** 2
    Ts    = 1.0 - Rs

    theta_t_nonTIR = np.rad2deg(np.arcsin(sin_theta_t[non_TIR]))
    return Rp, Tp, Rs, Ts, n1, n2, theta_t_nonTIR


# ── base class ────────────────────────────────────────────────────────────────

class MeasurementBase:

    def _init_arrays(self, tm_wavenum, te_wavenum, tm_single_beam, te_single_beam, nuextrema):
        """Set all shared spectral attributes. Called by subclass __init__ after any
        polarization-specific processing (e.g. Fresnel correction) is done."""
        self.TM_wavenum     = tm_wavenum
        self.TE_wavenum     = te_wavenum
        self.TM_single_beam = tm_single_beam
        self.TE_single_beam = te_single_beam

        self.TM_masked         = None
        self.TE_masked         = None
        self.TM_wavenum_masked = None
        self.TE_wavenum_masked = None

        if nuextrema is not None:
            mask = (te_wavenum > nuextrema[0]) & (te_wavenum < nuextrema[1])
            self.TE_masked         = te_single_beam[mask]
            self.TM_masked         = tm_single_beam[mask]
            self.TM_wavenum_masked = tm_wavenum[mask]
            self.TE_wavenum_masked = te_wavenum[mask]

    def alpha_ISB(self, background):
        """Return (wavenum, alpha_ISB) computed relative to a background measurement."""
        offset = np.log(background.TM_masked / background.TE_masked)
        alpha  = -np.log(self.TM_masked / self.TE_masked) + offset
        return self.TE_wavenum_masked, alpha

    def fit_lorentz(self, nu_range, nu0_guess, kappa_guess, background,
                    A_guess=None, B_guess=None):
        wavenum, alpha = self.alpha_ISB(background)
        return fit_lorentz(nu_range, nu0_guess, kappa_guess, wavenum, alpha,
                           A_guess, B_guess)

    def fit_nlorentz(self, nu_range, nu0_guesses, kappa_guesses, background,
                     A_guess=None, B_guesses=None):
        wavenum, alpha = self.alpha_ISB(background)
        return fit_nlorentz(nu_range, nu0_guesses, kappa_guesses, wavenum, alpha,
                            A_guess, B_guesses)

    def plot_single_beam(self, axs):
        """Plot raw TM and TE single-beam spectra."""
        axs.plot(self.TE_wavenum, self.TE_single_beam, label=f'TE {self.name}')
        axs.plot(self.TM_wavenum, self.TM_single_beam, label=f'TM {self.name}')

    def plot_ratio(self, axs):
        """Plot TM/TE polarization ratio over the masked wavenumber range."""
        axs.plot(self.TE_wavenum_masked, self.TM_masked / self.TE_masked,
                 label=f'TM/TE {self.name}')

    def plot_alpha_ISB(self, background, axs, label=None,add_label=True,lw=1.0):
        """Compute and plot alpha_ISB relative to background."""
        wavenum, alpha = self.alpha_ISB(background)
        if label is None:
            label = (
                rf"$-\ln \left(\frac{{I_{{{self.name},TM}}}}{{I_{{{self.name},TE}}}}\right)"
                rf"+ \ln \left(\frac{{I_{{{background.name},TM}}}}{{I_{{{background.name},TE}}}}\right)$"
            )
        if add_label:
            axs.plot(wavenum, alpha, label=label, linewidth=lw)
        else:
            axs.plot(wavenum, alpha, linewidth=lw)

    def plot_fit(self, fit_result, axs, nu_fit_plot_range=None, lw=1.0,
                component_colors=None, combined_color=None):
        """Plot a FitResult over this measurement's masked wavenumber range."""
        plot_lorentz_fit(fit_result, self.TE_wavenum_masked, axs, nu_fit_plot_range,
                         lw=lw, component_colors=component_colors, combined_color=combined_color)


# ── measurement subclasses ────────────────────────────────────────────────────

class MultipassMeas(MeasurementBase):
    def __init__(self, TEfile, TMfile, name, nuextrema=None):
        _, tm_wavenum, tm_single_beam, _ = load_data(TMfile)
        _, te_wavenum, te_single_beam, _ = load_data(TEfile)
        self.name = name
        self._init_arrays(tm_wavenum, te_wavenum, tm_single_beam, te_single_beam, nuextrema)


class SinglePassMeas(MeasurementBase):
    def __init__(self, TEfile, TMfile, name, thetai,
                 fresnel=False, n1=None, n2=None, n3=None, nuextrema=None):
        _, tm_wavenum, tm_single_beam, _ = load_data(TMfile)
        _, te_wavenum, te_single_beam, _ = load_data(TEfile)

        self.name              = name
        self.thetai            = thetai
        self.TM_single_beam_raw = tm_single_beam
        self.TE_single_beam_raw = te_single_beam

        if fresnel:
            _, Tp12, _, Ts12, _, _, theta_t12 = calculate_Fresnels(thetai, n1, n2)
            print("TP12:", Tp12, "Ts12:", Ts12)
            _, Tp23, _, Ts23, _, _, _ = calculate_Fresnels(theta_t12, n2, n3)
            print("TP23:", Tp23, "Ts23:", Ts23)
            tm_single_beam = tm_single_beam / (Tp12 * Tp23)
            te_single_beam = te_single_beam / (Ts12 * Ts23)

        self._init_arrays(tm_wavenum, te_wavenum, tm_single_beam, te_single_beam, nuextrema)


# ── backwards-compatible wrappers ─────────────────────────────────────────────

def fitLorentzPlot(nu_range, kappanu_guess, wavenum, alpha_ISB, axs_fits,
                   nu_fit_plot_range=None):
    alpha_select, wavenum_fit = maskFit(nu_range, wavenum, alpha_ISB)
    nu0_guess = wavenum_fit[np.argmax(alpha_select)]
    result    = fit_lorentz(nu_range, nu0_guess, kappanu_guess, wavenum, alpha_ISB)
    plot_lorentz_fit(result, wavenum, axs_fits, nu_fit_plot_range)
    return np.array([result.nu0s[0], result.kappas[0], result.A, result.Bs[0]])


def fitNormalPlot(nu_range, mu_guess, sigma_guess, wavenum, alpha_ISB, axs_fits):
    alpha_select, wavenum_fit = maskFit(nu_range, wavenum, alpha_ISB)

    A_guess   = alpha_select[-1]
    B_guess   = np.max(alpha_select) - np.min(alpha_select)
    fitGuess  = (mu_guess, sigma_guess, A_guess, B_guess)
    nu_width  = np.max(wavenum_fit) - np.min(wavenum_fit)
    fitBounds = ([np.min(wavenum_fit), 0,        -np.inf, 0      ],
                 [np.max(wavenum_fit), nu_width,  np.inf, np.inf ])

    fitNormal, _ = opt.curve_fit(fitFnNormal, wavenum_fit, alpha_select,
                                 p0=fitGuess, bounds=fitBounds)
    fit_label = (r"$\mu = %0.2f\ \mathrm{cm}^{-1},\ \sigma = %0.2f\ \mathrm{cm}^{-1}$"
                 % (fitNormal[0], fitNormal[1]))
    axs_fits.plot(wavenum, [fitFnNormal(nu, *fitNormal) for nu in wavenum],
                  linewidth=1, label=fit_label)
    return fitNormal
