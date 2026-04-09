import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
from FTIR_analysis_helpers import build_MP, fit_lorentz, fit_nlorentz, plot_lorentz_fit

# ── sample / settings ────────────────────────────────────────────────────────
sample_name          = 'P12-9-25-2MP'
background_samp_name = 'P12-9-25-1-GaAs-intrinsic-only-MP'
preamp               = 'PA101'

numin, numax = 825, 3100

base_dir = r'C:\Users\sp6497_a\OneDrive - Princeton University\20260205_P12-9-25-2-MP\20260211'
bg_dir   = r'C:\Users\sp6497_a\OneDrive - Princeton University\P12-9-25-1-intrinsic-GaAs-only-MP'

# ── load data ─────────────────────────────────────────────────────────────────
tm_file    = os.path.join(base_dir, sample_name + '-P0deg_novislight.CSV')
te_file    = os.path.join(base_dir, sample_name + '-P90deg_novislight.CSV')
tm_bg_file = os.path.join(bg_dir,   background_samp_name + '-P0deg.CSV')
te_bg_file = os.path.join(bg_dir,   background_samp_name + '-P90deg.CSV')

samp_meas = build_MP(te_file,    tm_file,    sample_name,          nuextrema=[numin, numax])
bg_meas   = build_MP(te_bg_file, tm_bg_file, background_samp_name, nuextrema=[numin, numax])

# ── ISB absorption ────────────────────────────────────────────────────────────
offset    = np.log(bg_meas.TM_masked / bg_meas.TE_masked)
alpha_ISB = -np.log(samp_meas.TM_masked / samp_meas.TE_masked) + offset

wavenum = samp_meas.TE_wavenum_masked

# ── raw scans + ratios figure ─────────────────────────────────────────────────
fig, axs = plt.subplots(1, 3, figsize=(14, 8))
fig.suptitle(sample_name)

axs[0].set_title('Raw data')
axs[0].set_xlabel('Wavenumber (cm^-1)', fontsize=12)
axs[0].set_ylabel('Single Beam', fontsize=12)

snr_title = f'SNR mask {numin}' + r'$ < \nu < $' + f'{numax}' + r' ${cm}^{-1}$'
axs[1].set_title(snr_title)
axs[1].set_xlabel('Wavenumber (cm^-1)', fontsize=12)
axs[1].set_ylabel('Polarization Transmission Ratio', fontsize=12)

axs[2].set_title(snr_title)
axs[2].set_xlabel('Wavenumber (cm^-1)', fontsize=12)
axs[2].set_ylabel('Samp over bg transmission ratios', fontsize=12)

for ax in axs:
    ax.grid()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

axs[0].plot(samp_meas.TE_wavenum, samp_meas.TE_single_beam, color='blue',  label='TE ' + samp_meas.name)
axs[0].plot(samp_meas.TM_wavenum, samp_meas.TM_single_beam, color='red',   label='TM ' + samp_meas.name)
axs[0].plot(bg_meas.TE_wavenum,   bg_meas.TE_single_beam,   color='c',     label='TE ' + bg_meas.name)
axs[0].plot(bg_meas.TM_wavenum,   bg_meas.TM_single_beam,   color='m',     label='TM ' + bg_meas.name)

axs[1].plot(samp_meas.TM_wavenum_masked, samp_meas.TM_masked / samp_meas.TE_masked, label='TM/TE ' + sample_name)
axs[1].plot(bg_meas.TM_wavenum_masked,   bg_meas.TM_masked   / bg_meas.TE_masked,   label='TM/TE ' + bg_meas.name)

axs[2].plot(samp_meas.TM_wavenum_masked, samp_meas.TM_masked / bg_meas.TM_masked, label='TM samp / TM bg')
axs[2].plot(samp_meas.TE_wavenum_masked, samp_meas.TE_masked / bg_meas.TE_masked, label='TE samp / TE bg')

for ax in axs:
    ax.legend(prop={'size': 10})
plt.tight_layout()
plt.savefig(os.path.join(base_dir, sample_name + '_raw_scans_and_ratios.svg'))

# ── alpha_ISB + fits figure ───────────────────────────────────────────────────
fig_fits, axs_fits = plt.subplots(figsize=(10, 8))
axs_fits.set_title(sample_name + f' SNR mask {numin}' + r'$ < \nu < $' + f'{numax}')
axs_fits.set_xlabel('Wavenumber (cm^-1)', fontsize=12)
axs_fits.set_ylabel(r'$\alpha_{ISB} \times L_{path}$', fontsize=12)
axs_fits.grid()
axs_fits.xaxis.set_major_locator(MaxNLocator(integer=True))

axs_fits.plot(
    wavenum, alpha_ISB, color='black',
    label=rf"$-\ln \left(\frac{{I_{{{samp_meas.name},TM}}}}{{I_{{{samp_meas.name},TE}}}}\right)"
          rf"+ \ln \left(\frac{{I_{{{bg_meas.name},TM}}}}{{I_{{{bg_meas.name},TE}}}}\right)$"
)

# ── Option A: fit each peak independently ─────────────────────────────────────
# result_peak1 = fit_lorentz(
#     nu_range    = [948, 1211],
#     nu0_guess   = 1080,   # adjust to your expected peak position
#     kappa_guess = 50,
#     wavenum     = wavenum,
#     alpha_ISB   = alpha_ISB,
# )
# plot_lorentz_fit(result_peak1, wavenum, axs_fits)

# result_peak2 = fit_lorentz(
#     nu_range    = [1270, 1594],
#     nu0_guess   = 1430,   # adjust to your expected peak position
#     kappa_guess = 100,
#     wavenum     = wavenum,
#     alpha_ISB   = alpha_ISB,
# )
# plot_lorentz_fit(result_peak2, wavenum, axs_fits)

# ── Option B: fit both peaks simultaneously (comment out A above to compare) ──
result_both = fit_nlorentz(
    nu_range      = [948, 1594],
    nu0_guesses   = [1086, 1300],
    kappa_guesses = [50, 100],
    wavenum       = wavenum,
    alpha_ISB     = alpha_ISB,
)
plot_lorentz_fit(result_both, wavenum, axs_fits)

axs_fits.legend(prop={'size': 14})
plt.tight_layout()
plt.savefig(os.path.join(base_dir, sample_name + '_alpha_ISB.svg'))

plt.show()
