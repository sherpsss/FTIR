
import pandas as pd
import matplotlib.pyplot as plt
import os

def load_data(filename):
    data = pd.read_csv(filename, header=None)
    wavenumber = data[0]  # First column: wavenumber
    single_beam = data[1]  # Second column: transmission
    wavelength = 1e4 / wavenumber  # Convert wavenumber (cm^-1) to wavelength (µm)

    # Extract a label from the filename (remove path and extension)
    label = os.path.splitext(os.path.basename(filename))[0]

    return wavelength, wavenumber, single_beam, label

class ThetaIncident:
    def __init__(self,angle):
        self.angle = angle
        self.TM_single_beam = None
        self.TE_single_beam = None
        self.TM_wavenum = None
        self.TE_wavenum = None

        self.TE_background_single_beam = None
        self.TM_background_single_beam = None
        self.TM_bg_wavenum = None
        self.TE_bg_wavenum = None

        self.TM_ratio = None
        self.TE_ratio = None

        
def build_theta_dict(base_dir,tm_suffix,te_suffix):
    
    theta_dict = {}

    background_te_file = os.path.join(base_dir,'backgrounds/background_'+te_suffix+'.CSV')
    background_tm_file = os.path.join(base_dir,'backgrounds/background_'+tm_suffix+'.CSV')

    _,background_tm_wavenum, background_tm_single_beam,_ = load_data(background_tm_file)
    _,background_te_wavenum,background_te_single_beam,_ = load_data(background_te_file)

    for folder_name in os.listdir(base_dir):

        if folder_name.startswith('theta_i_'):
            theta_key_deg_str = folder_name.split('_')[-1]
            theta_key = theta_key_deg_str.split('deg')[0]
            theta_path = os.path.join(base_dir,folder_name)

            theta_i_obj = ThetaIncident(angle = theta_key)
            theta_i_obj.TM_bg_wavenum = background_tm_wavenum
            theta_i_obj.TE_bg_wavenum = background_te_wavenum
            theta_i_obj.TM_background_single_beam = background_tm_single_beam
            theta_i_obj.TE_background_single_beam = background_te_single_beam

            #now fill in the params for theta_incident (TM vs TE)
            for filename in os.listdir(theta_path):
                if tm_suffix in filename:
                    tm_theta_i_filepath = os.path.join(theta_path,filename)
                    _, tm_wavenum, tm_single_beam, _ = load_data(tm_theta_i_filepath)
                    theta_i_obj.TM_single_beam = tm_single_beam
                    theta_i_obj.TM_ratio = tm_single_beam/background_tm_single_beam
                    theta_i_obj.TM_wavenum = tm_wavenum
                elif te_suffix in filename:
                    te_theta_i_filepath = os.path.join(theta_path,filename)
                    _, te_wavenum, te_single_beam, _ = load_data(te_theta_i_filepath)
                    theta_i_obj.TE_single_beam = te_single_beam
                    theta_i_obj.TE_ratio = te_single_beam/background_te_single_beam
                    theta_i_obj.TE_wavenum = te_wavenum


            theta_dict[theta_key]=theta_i_obj
    return theta_dict
            
te_polarization = 'P90deg'
tm_polarization = 'P0deg'

samp_dir = '/Users/srsplatt/Library/Mobile Documents/com~apple~CloudDocs/Princeton/Gmachl Research/20241211_theta_i_sweep_singlepass_342A'

theta_dict = build_theta_dict(samp_dir,tm_polarization,te_polarization)

fig, axs = plt.subplots(2, 1, figsize=(10, 10))
axs[0].set_xlabel('Wavenumber (cm^-1)')
axs[0].set_ylabel(r"$T_{\theta_i}= I_{out,sample}/I_{out,background}$")
theta_variation_title = 'theta_i single pass 342A'
axs[0].set_title(theta_variation_title)

axs[1].set_xlabel('Wavenumber (cm^-1)')
axs[1].set_ylabel(r"$T_{\theta_i,TM}/T_{\theta_i,TE}$")
theta_variation_title_polarization_ratios = theta_variation_title + ' polarization ratios'
axs[1].set_title(theta_variation_title_polarization_ratios)
#iterate through the dictionary

#TM TE ratios
# fig_Pratios, axs_Pratios = plt.subplots(1,1,figsize=(12,8))
# axs_Pratios.set_xlabel('Wavenumber (cm^-1)')
# axs_Pratios.set_ylabel('polarization transmission ratios [TE/TM]')
# polarization_variation_titles = theta_variation_title + 'TE TM ratios'
# axs_Pratios.set_title(polarization_variation_titles)

colors_dict = {}
colors_dict['0'] = 'blue'
colors_dict['15'] = 'green'
colors_dict['30'] = 'red'
colors_dict['40'] = 'orange'

for angle in theta_dict.keys():
    theta_i_data = theta_dict[angle]

    #generate the specific angle plot for sanity check
    fig_single_theta, ax_single_theta = plt.subplots(1, 1, figsize=(12, 8))
    ax_single_theta.scatter(theta_i_data.TE_wavenum,theta_i_data.TE_single_beam,label = r"$\theta_i$ = " + angle + 'TE',color='blue',marker='o',s=7)
    ax_single_theta.scatter(theta_i_data.TM_wavenum,theta_i_data.TM_single_beam,label = r"$\theta_i$ = " + angle + 'TM',color='red',marker='o',s=7)
    ax_single_theta.scatter(theta_i_data.TM_bg_wavenum,theta_i_data.TM_background_single_beam,label = 'background TM',color='red',marker='x',s=15)
    ax_single_theta.scatter(theta_i_data.TE_bg_wavenum,theta_i_data.TE_background_single_beam,label='background TE',color='blue',marker='x',s=15)

    ax_single_theta.set_title('raw data' + angle + r"= $\theta_i$")
    ax_single_theta.set_xlabel(r"wavenumber [${cm}^{-1}$]")
    ax_single_theta.set_ylabel("Single Beam")
    ax_single_theta.legend()
    ratio_label = f'{angle} raw data'
    save_title = ratio_label + '.svg'
    # plt.show()
    fig_single_theta.savefig(os.path.join(samp_dir, save_title))
    # fig_single_theta.show()
    #generate additions for the summary plots
    axs[0].scatter(theta_i_data.TE_wavenum,theta_i_data.TE_ratio,label = r"$\theta_i$ = " + angle + 'TE',color=colors_dict[angle],marker ='o',s=2)
    axs[0].scatter(theta_i_data.TM_wavenum,theta_i_data.TM_ratio,label = r"$\theta_i$ = " + angle + 'TM',color=colors_dict[angle],marker='x',s=4)

    #now consider the TM TE ratios
    TE_TM_ratio = theta_i_data.TM_ratio/theta_i_data.TE_ratio
    axs[1].scatter(theta_i_data.TE_wavenum,TE_TM_ratio,label = r"$\theta_i$ = " + angle,color=colors_dict[angle],s=2)
    # axs.scatter(theta_i_data.TM_wavenum,theta_i_data.TM_ratio,label = r"$\theta_i$ = " + angle + 'TM',color=colors_dict[angle],marker='x',s=15)


save_title_comparison_plot = theta_variation_title + '.svg'
axs[0].legend()
axs[1].legend()
plt.figure(fig)
plt.show()
fig.savefig(os.path.join(samp_dir, save_title_comparison_plot))
