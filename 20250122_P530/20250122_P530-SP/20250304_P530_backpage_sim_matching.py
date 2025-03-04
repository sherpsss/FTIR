import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

def plot_transitions(stateindices_arr,deltanus_arr,dipoles_arr,botstatenum_arr,base_dir,wafername,expected_trans=None,botstatenum_markers=None,titlesize=20,ticksize=12,axislabelsfont=12,markersize=10):
    fig, axs = plt.subplots(figsize=(10,8))
    axs_dipole = axs.twinx()
    axs.set_xlabel('state n', fontsize=axislabelsfont)
    axs.grid(True)
    axs_dipole.set_ylabel(r"$|d_{ij}|[A]$", fontsize=axislabelsfont,color='b')
    axs.set_ylabel(r"$\Delta E[meV]$",fontsize=axislabelsfont,color='r')
    marker = "."
    # axs.set_title(wafername + " simulated ISB transitions from state n to state " + str(botstatenum))
    for botstatenum_ind in range(0,len(botstatenum_arr)):
        stateindices = stateindices_arr[botstatenum_ind]
        deltanus = deltanus_arr[botstatenum_ind]
        botstatenum = botstatenum_arr[botstatenum_ind]
        dipoles = dipoles_arr[botstatenum_ind]
        if botstatenum_markers is not None:
            marker = botstatenum_markers[botstatenum_ind]
        axs.scatter(stateindices, deltanus, label=r"$\nu$"+str(botstatenum)+r"$-\nu n$",color='r',marker=marker,s=markersize)
        axs_dipole.scatter(stateindices, np.abs(dipoles), label=r"$d$" + str(botstatenum) + r"$n$",color='b',marker=marker,s=markersize)


    if expected_trans is not None:
        measured_nu = np.ones([len(stateindices_arr[0])])*expected_trans[0]
        lower_range = np.ones([len(stateindices_arr[0])])*(expected_trans[1]/2)
        axs.plot(stateindices_arr[0],measured_nu, color='g',label=r"${\Delta E}_{measured}$",alpha=0.6)
        axs.fill_between(stateindices_arr[0],measured_nu-lower_range,measured_nu+lower_range,color='g',label=r"$\Gamma_{measured}$",alpha=0.3)
        axs.fill_between(stateindices_arr[0],measured_nu*0.9,measured_nu*1.1,color='royalblue',label="sim. tol.",alpha=0.3)
        # axs.axhline()
    # axs.legend()
    # axs.legend(prop={"size": 14})
    # axs_dipole.legend()
    # plt.tight_layout()
    # save_title = os.path.join(base_dir, wafername + ' E' + str(botstatenum) + 'to En' + '.svg')

    axs.tick_params(axis='x', labelsize=ticksize)
    axs.tick_params(axis='y', labelsize=ticksize)
    axs_dipole.tick_params(axis='y', labelsize=ticksize)
    save_title = os.path.join(base_dir, wafername + ' Ecombined' + '.svg')
    plt.savefig(save_title)
    return fig,axs,axs_dipole
meV_per_inv_cm = 8.065
wafername = 'P530'
base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250122_P530-A-SP/simulations'
simfilename = 'P530 transitions - consistent exports - P530 backpage'
sim_file = os.path.join(base_dir, simfilename + '.CSV')

data = pd.read_csv(sim_file, header=None)
statenum_repetitive = np.array(data[0])
statenum = np.array(data[1])  # First column: statenum
deltaE0En_meV = np.array(data[2])  # column 1: transition gap state 2->n
deltaE1En_meV = np.array(data[4]) # column 3: transition gap state 5->n

dipole0n = np.array(data[3]) # column 2: dipole transition gap state 2 -> n
dipole1n = np.array(data[5]) #column 4: dipole transition gap state 5-> n

# deltanu0n = deltaE0En_meV*meV_per_inv_cm
# deltanu1n = deltaE1En_meV*meV_per_inv_cm

blues = ['darkblue','mediumblue','blue','cornflowerblue']
reds = ['mediumvioletred','deeppink','hotpink','pink']

measured_deltanu = 1127.0
measured_deltanu_meV = measured_deltanu/meV_per_inv_cm
measured_gamma = 92.23
measured_gamma_meV = measured_gamma/meV_per_inv_cm

stateindices_arr = [statenum[1:],statenum[1:]]
# deltanus_arr = [deltanu0n[1:],deltanu1n[1:]]
deltanus_arr = [deltaE0En_meV[1:],deltaE1En_meV[1:]]
dipoles_arr = [dipole0n[1:],dipole1n[1:]]
botstatenum_arr = [1,2]

axislabelsfont=30
legendfont=30
ticksize=30
titlesize=30
markersize = 300
# gnd_state_plot,gnd_state_ax = plot_transitions(statenum,deltanu0n,dipole0n,0,base_dir,wafername,expected_trans=[measured_deltanu,measured_gamma])
gnd_state_plot,gnd_state_ax,gnd_state_dipole_ax = plot_transitions(stateindices_arr,deltanus_arr,dipoles_arr,botstatenum_arr,base_dir,wafername,expected_trans=[measured_deltanu_meV,measured_gamma_meV],botstatenum_markers=["+","."],titlesize=titlesize,ticksize=ticksize,axislabelsfont=axislabelsfont,markersize=markersize)
# E1state_plot,E1_state_ax = plot_transitions(statenum,deltanu1n,dipole1n,1,base_dir,wafername,expected_trans=[measured_deltanu,measured_gamma])

handles, labels = gnd_state_ax.get_legend_handles_labels()
handles_dipole, labels_dipole = gnd_state_dipole_ax.get_legend_handles_labels()
by_label = dict(zip(labels+labels_dipole, handles+handles_dipole))
# plt.legend(by_label.values(), by_label.keys())
gnd_state_ax.legend(by_label.values(), by_label.keys(),prop={"size": legendfont})

plt.figure(gnd_state_plot)
plt.tight_layout()
save_title = os.path.join(base_dir, wafername + ' Ecombined_legend' + '.svg')
plt.savefig(save_title)

plt.figure(gnd_state_plot)
plt.show()


# plt.figure(E1state_plot)
# plt.show()