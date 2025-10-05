import numpy as np
import matplotlib.pyplot as plt

n_air = 1.0
n_InP = 2.7132  # InP at lambda = 8 um
n_GaAs = 3.28
n_ZnSe = 2.39


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



def plot_Fresnel(theta_is, Rp, Tp, Rs, Ts, n1, n2,theta_ts):

    fig, axs = plt.subplots(1, 2, figsize=(14, 8))
    axs[0].plot(theta_is, Rp, label=r"$R_p$", color="blue")
    axs[0].plot(theta_is, Tp, label=r"$T_p$", color="blue", linestyle=":")
    axs[0].plot(theta_is, Rs, label=r"$R_s$", color="red")
    axs[0].plot(theta_is, Ts, label=r"$T_s$", color="red", linestyle=":")

    axs[0].set_xlabel(r"Incident Angle $\theta_i$ (degrees)")
    axs[0].set_ylabel(r"Power coefficients")
    axs[0].set_title(r"power coefficients vs. $\theta_i$, $n_1$ = " + str(n1) + ", $n_2$ = " + str(n2))
    axs[0].grid(True)
    axs[0].legend()

    if theta_ts is not None:

        axs[1].plot(theta_is[0:len(theta_ts)],theta_ts,label=r"$theta_t")
        axs[1].set_xlabel(r"Incident angle $\theta_i$ (degrees)")
        axs[1].set_ylabel(r"Transmission angle $\theta_i$ (degrees)")

    return fig,axs


theta_is = np.linspace(0, 90, 500)
# first interface (coupling in)
# Rp_A, Tp_A, Rs_A, Ts_A, n1_A, n2_A = calculate_Fresnels(theta_is, n_air, n_InP)
#
# plot_Fresnel(theta_is, Rp_A, Tp_A, Rs_A, Ts_A, n1_A, n2_A)
# second interface (InP to air polished facet)

Rp_B, Tp_B, Rs_B, Ts_B,n1_B, n2_B,theta_ts = calculate_Fresnels(theta_is, n_air, n_GaAs)
fig_ent,axs_ent=plot_Fresnel(theta_is, Rp_B, Tp_B, Rs_B, Ts_B, n1_B, n2_B,theta_ts)

Rp_B, Tp_B, Rs_B, Ts_B, n1_B, n2_B,theta_ts = calculate_Fresnels(theta_is, n_GaAs, n_air)
fig_ex,axs_ex=plot_Fresnel(theta_is, Rp_B, Tp_B, Rs_B, Ts_B, n1_B, n2_B,theta_ts)

Rp_B, Tp_B, Rs_B, Ts_B, n1_B, n2_B,theta_ts = calculate_Fresnels(theta_is, n_ZnSe, n_GaAs)
fig_ATR,axs_ATR=plot_Fresnel(theta_is, Rp_B, Tp_B, Rs_B, Ts_B, n1_B, n2_B,theta_ts)

plt.figure(fig_ent)
plt.show()

plt.figure(fig_ex)
plt.show()

plt.figure(fig_ATR)
plt.show()