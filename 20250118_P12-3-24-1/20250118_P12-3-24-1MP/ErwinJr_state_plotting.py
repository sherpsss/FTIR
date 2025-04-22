import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

class QuantizedState:
    def __init__(self,num):
        self.Energy = None
        self.WF_ys = None
        self.WF_xs = None

def generate_MQW_plot(base_dir,base_filename,title,states_to_plot=None,axislabelsfont=18,legendfont=18,ticksize=18,titlesize=20):
    wavefunctions_file = os.path.join(base_dir, base_filename+'_WFs.csv')
    conduction_band_file = os.path.join(base_dir, base_filename+'_CB.csv')
    energies_file = os.path.join(base_dir, base_filename+'_Es.csv')

    # Load wavefunction data and ensure numeric values
    wavefunctions_str = pd.read_csv(wavefunctions_file, header=None)
    wf_data = wavefunctions_str.apply(pd.to_numeric)  # Convert all values to floats
    wf_xs = np.array(wf_data[0])

    conduction_bands_str = pd.read_csv(conduction_band_file, header=None)
    cb_data = conduction_bands_str.apply(pd.to_numeric)  # Convert all values to floats
    cb_xs = np.array(cb_data[0])
    cb_ys = np.array(cb_data[1])
    print("minimum conduction band value: " + str(cb_ys[0]))

    Es_str = pd.read_csv(energies_file, header=None)
    Es_data = Es_str.apply(pd.to_numeric)
    Es = np.array(Es_data[0])

    # Plot
    fig_Es, axs_Es = plt.subplots(figsize=(10, 10))
    # axs_Es.plot(cb_xs, cb_ys, label="conduction band", color='black')
    axs_Es.plot(cb_xs, cb_ys, color='black')


    axs_Es.tick_params(axis='x', labelsize=ticksize)
    axs_Es.tick_params(axis='y', labelsize=ticksize)

    axs_Es.set_xlabel(r"growth direction [Angstroms]", fontsize=axislabelsfont)
    axs_Es.set_ylabel(r"Energy offset from ${Ga}_{7}{As}_{10}$ CBM[eV]", fontsize=axislabelsfont)

    axs_Es.set_title(title,fontsize=titlesize)

    # build the state dict
    # states_to_plot = np.array([0, 1, 5, 6, 10, 11])
    if states_to_plot is None:
        states_to_plot = np.arange(len(Es))

    for E_ind in states_to_plot:
        E_state = QuantizedState(num=E_ind)
        E_state.WF_xs = wf_xs
        E_state.WF_ys = np.array(wf_data[E_ind + 1])
        E_state.Energy = Es[E_ind]
        if E_ind ==1 or E_ind==3:
            print("state " + str(E_ind) + " energy = " + str(E_state.Energy))


        #     wavefunctions_normalized.iloc[:, i] = wavefunctions.iloc[:, i] / np.max(np.abs(wavefunctions.iloc[:, i]))
        #     wavefunctions_normalized.iloc[:, i] *= 0.1  # Scale factor for visibility
        # axs_Es.plot(E_state.WF_xs, (E_state.WF_ys) * 0.3 + E_state.Energy, label=str(E_state.Energy))
        axs_Es.plot(E_state.WF_xs, (E_state.WF_ys) * 0.3 + E_state.Energy)
    save_title = os.path.join(base_dir, base_filename + '.svg')

    # plt.legend()
    # axs_Es.legend(prop={"size": legendfont})
    plt.tight_layout()
    plt.savefig(save_title)
    return axs_Es,fig_Es



# Load data
GaAs_width_nm_larger = 9
GaAs_width_nm_smaller = 7
sim_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/P12-3-24-1/simulations'
sim_dir_P530 = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250122_P530-A-SP/simulations'

GaAs_larger_filename = "P12-3-24-1_stack_933nu_transition"
GaAs_smaller_filename = "P12-3-24-1_stack_7nm_GaAs_well"
P530_filename = "P530_stack_2repeats"

states_to_plot_GaAs_larger = np.array([0,1,5,6,10,11])
states_to_plot_GaAs_smaller = np.array([0,1,5,6,10,11])
# states_to_plot_GaAs_smaller = np.array([1,6,11])
# states_to_plot_P530 = np.array([1,3,12,14])

axislabelsfont=20
legendfont=20
ticksize=20
titlesize=20

# axes_GaAs_larger,fig_GaAs_larger = generate_MQW_plot(sim_dir,GaAs_larger_filename,"Example Sample, GaAs well = " + str(GaAs_width_nm_larger) + " nm, 2 periods",states_to_plot=states_to_plot_GaAs_larger,axislabelsfont=axislabelsfont,legendfont=legendfont,ticksize=ticksize,titlesize=titlesize)
axes_GaAs_smaller,fig_GaAs_smaller = generate_MQW_plot(sim_dir,GaAs_smaller_filename,"P2-5-25-1, 2 periods",states_to_plot=states_to_plot_GaAs_smaller,axislabelsfont=axislabelsfont,legendfont=legendfont,ticksize=ticksize,titlesize=titlesize)

# axes_P530,fig_P530 = generate_MQW_plot(sim_dir_P530,P530_filename,"Complex Sample, 2 periods",states_to_plot=states_to_plot_P530,axislabelsfont=axislabelsfont,legendfont=legendfont,ticksize=ticksize,titlesize=titlesize)
# axes_P530,fig_P530 = generate_MQW_plot(sim_dir_P530,P530_filename,"Complex Sample, 2 periods",axislabelsfont=axislabelsfont,legendfont=legendfont,ticksize=ticksize,titlesize=titlesize)
# plt.figure(fig_P530)
# plt.tight_layout()
# plt.show()
#
# axes_GaAs_larger.set_xlim(xmin=1100,xmax=1700)
period_smaller = 320
target_smaller = 3

axes_GaAs_smaller.set_xlim(xmin=period_smaller*3,xmax=period_smaller*4+period_smaller)

# plt.figure(fig_GaAs_larger)
# plt.tight_layout()
# plt.savefig(os.path.join(sim_dir, GaAs_larger_filename + 'xlimed.svg'))

plt.figure(fig_GaAs_smaller)
plt.tight_layout()
plt.savefig(os.path.join(sim_dir, GaAs_smaller_filename + 'xlimed.svg'))


# plt.figure(fig_GaAs_larger)
# plt.show()
plt.figure(fig_GaAs_smaller)
plt.show()
