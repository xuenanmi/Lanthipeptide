import glob
import matplotlib.pyplot as plt
import numpy as np
import pickle
from matplotlib.ticker import FormatStrFormatter

WT_totdist = []
for file in sorted(glob.glob('/home/xmi4/RippsProject/procA33_WT/analysis/*all-dist-per-traj.npy')):
    distI = np.load(file)
    WT_totdist.append(distI)
WT_totdist = np.concatenate(WT_totdist)
WT_weights = pickle.load(open('/home/xmi4/RippsProject/procA33_WT/analysis/MSM-procA3.3WT_cluster_150_ticdim_10-weights.pkl','rb'))
WT_weights = np.concatenate(WT_weights)

MT1_totdist = []
for file in sorted(glob.glob('/home/xmi4/RippsProject/procA33_MT1_2/analysis/*all-dist-per-traj.npy')):
    distI = np.load(file)
    MT1_totdist.append(distI)
MT1_totdist = np.concatenate(MT1_totdist)
MT1_weights = pickle.load(open('/home/xmi4/RippsProject/procA33_MT1_2/analysis/MSM-procA3.3MT1_cluster_150_ticdim_8-weights.pkl','rb'))
MT1_weights = np.concatenate(MT1_weights)


MT2_totdist = []
for file in sorted(glob.glob('/home/xmi4/RippsProject/procA33_MT2_16/analysis/*all-dist-per-traj.npy')):
    distI = np.load(file)
    MT2_totdist.append(distI)
MT2_totdist = np.concatenate(MT2_totdist)
MT2_weights = pickle.load(open('/home/xmi4/RippsProject/procA33_MT2_16/analysis/MSM-procA3.3MT2_cluster_100_ticdim_6-weights.pkl','rb'))
MT2_weights = np.concatenate(MT2_weights)

def plot_hist(feat, bins_, title_, xlabel_, xticks_, yticks_):
    nSD, binsSD = np.histogram(WT_totdist[:,feat]*10, bins=bins_, density=True, weights = WT_weights)
    nSD1, binsSD1 = np.histogram(MT1_totdist[:,feat]*10, bins=bins_, density=True, weights = MT1_weights)
    nSD2, binsSD2 = np.histogram(MT2_totdist[:,feat]*10, bins=bins_, density=True, weights = MT2_weights)
    #averageSD = [(binsSD[j]+binsSD[j+1])/2 for j in range(len(binsSD)-1)]

    fig, axs = plt.subplots(1,1,figsize=(10,7))
    axs.plot(binsSD[0:-1], nSD, linewidth=3, c='red')
    axs.plot(binsSD1[0:-1], nSD1, linewidth=3, c='blue')
    axs.plot(binsSD2[0:-1], nSD2, linewidth=3, c='orange')
    axs.legend(['Wild Type', 'Mutant1', 'Mutant2'], fontsize = 18)
    
    #axs.set_xticks(np.arange(int(min(binsSD[0:-1])),int(max(binsSD[0:-1]))+1,1))
    #axs.set_xticklabels(np.arange(int(min(binsSD[0:-1])),int(max(binsSD[0:-1]))+1,1))
    axs.set_xticks(xticks_)
    axs.set_xticklabels(xticks_)

    axs.set_yticks(yticks_)
    axs.set_yticklabels(yticks_)
    axs.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axs.tick_params(width=3,length=5, labelsize=20)
    
    plt.xlabel(xlabel_, fontsize=28)
    plt.ylabel('Probability Density', fontsize=28)
    plt.rc('xtick', labelsize=24)
    plt.rc('ytick', labelsize=24)
    
    plt.tight_layout()
    plt.savefig('Prob_dist_plot_' + title_ + '.png',dpi=500)

#plot_hist(155, 100, '11_14')
#plot_hist(162, 100, '11_21')
#plot_hist(183, 100, '14_18')
#plot_hist(204, 100, '18_21')

#plot_hist(155, 100,'11_14', 'T11-C14 Distance ($\AA$)', np.arange(0,11,2), np.arange(0,0.81,0.1))
plot_hist(162, 100, '11_21', 'T11-C21 Distance ($\AA$)', np.arange(0,31,5), np.arange(0,0.351,0.05))
#plot_hist(183, 100, '14_18', 'C14-T18 Distance ($\AA$)', np.arange(0,15,2), np.arange(0,0.41,0.1))
#plot_hist(204, 100, '18_21', 'T18-C21 Distance ($\AA$)', np.arange(0,11,2), np.arange(0,0.81,0.1))

