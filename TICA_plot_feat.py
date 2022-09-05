import matplotlib.pyplot as plt
import mdtraj as md
import os
import glob
import numpy as np
import pickle
import pyemma
import seaborn as sns

totdist = []
for file in glob.glob('*all-dist-per-traj.npy'):
    distI = np.load(file)
    #distI = np.delete(distI,loop,axis=1)
    #print(file)
    totdist.append(distI)

data_tic = pyemma.coordinates.tica(totdist,lag=250, dim=10).get_output()
data_tic_concatenated = np.concatenate(data_tic)
totdist_concatenated = np.concatenate(totdist)

fig, axs = plt.subplots(1,1,figsize=(10,7))
sc = plt.scatter(data_tic_concatenated[:,0], data_tic_concatenated[:,1],c=totdist_concatenated[:,155]*10,vmin =np.min(totdist_concatenated[:,155]*10),vmax = np.max(totdist_concatenated[:,155]*10),cmap='jet')
cbar = fig.colorbar(sc,ticks=range(int(np.max(totdist_concatenated[:,155]*10))+2))
cbar.set_label("Thr11-Cys14 Distance ($\AA$)",size=24)
cbar.ax.tick_params(labelsize=16)


plt.xlabel('Projection on 1st tic', fontsize=24)
plt.ylabel('Projection on 2nd tic', fontsize=24)

plt.tight_layout()
plt.savefig('WT_tica_feat_T11_C14',dpi=500)
