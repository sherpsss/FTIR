import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

def plot_transitions(stateindices,deltanus,dipoles,botstatenum,base_dir,wafername):
    fig, axs = plt.subplots(figsize=(12, 8))
    axs.scatter(stateindices, deltanus, label=r"$\nu$"+str(botstatenum)+r"$-\nu n$",color='r')
    axs_dipole = axs.twinx()
    axs_dipole.scatter(stateindices, np.abs(dipoles), label=r"dipole strength " + str(botstatenum) + r" to n",color='b')

    axs.set_xlabel('state n', fontsize=12)
    axs.grid(True)
    axs_dipole.set_ylabel(r"$|d_{ij}|[A]$", fontsize=12,color='b')
    axs.set_ylabel(r"$\Delta \nu [{cm}^{-1}]$",fontsize=12,color='r')
    axs.set_title(wafername + " simulated ISB transitions from state n to state " + str(botstatenum))

    axs.legend()
    axs.legend(prop={"size": 14})
    axs_dipole.legend()
    plt.tight_layout()
    save_title = os.path.join(base_dir, wafername + ' E' + str(botstatenum) + 'to En' + '.svg')
    plt.savefig(save_title)
    return fig,axs

meV_per_inv_cm = 8.065
wafername = 'P530'
base_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20250122_P530-A-SP/simulations'
simfilename = 'P530 transitions - backpage growth'
sim_file = os.path.join(base_dir, simfilename + '.CSV')

data = pd.read_csv(sim_file, header=None)
statenum = np.array(data[0])  # First column: statenum
deltaE0En_meV = np.array(data[1])  # column 1: transition gap state 2->n
deltaE1En_meV = np.array(data[3]) # column 3: transition gap state 5->n

dipole0n = np.array(data[2]) # column 2: dipole transition gap state 2 -> n
dipole1n = np.array(data[4]) #column 4: dipole transition gap state 5-> n

deltanu0n = deltaE0En_meV*meV_per_inv_cm
deltanu1n = deltaE1En_meV*meV_per_inv_cm

blues = ['darkblue','mediumblue','blue','cornflowerblue']
reds = ['mediumvioletred','deeppink','hotpink','pink']

measured_deltanu = 1127.0

gnd_state_plot,gnd_state_ax = plot_transitions(statenum,deltanu0n,dipole0n,0,base_dir,wafername)
gnd_state_ax.axhline(y=measured_deltanu,color='g')

E1state_plot,E1_state_ax = plot_transitions(statenum,deltanu1n,dipole1n,1,base_dir,wafername)
E1_state_ax.axhline(y=measured_deltanu,color='g')

plt.figure(gnd_state_plot)
plt.show()

plt.figure(E1state_plot)
plt.show()