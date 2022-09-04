import numpy as np
import pandas as pd
import glob
import mdtraj as md
import math
import matplotlib.pyplot as plt
import os
from numpy import linalg as LA
import pickle
import pyemma
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc
import sys
hfont = {'fontname':'Helvetica'}


df = pd.DataFrame(columns = ['feat1','feat2','feat3','feat4'])
#df_MT1 = pd.DataFrame(columns = ['feat1','feat2','feat3','feat4'])
threshold = 0.45

"""
for j in range(200):
    files = pickle.load(open("/home/xmi4/RippsProject/procA1.1/plot/bootstrapping_ticdim/bt_80_" + str(j) +"_files.pkl",'rb'))
    msm = pickle.load(open("/home/xmi4/RippsProject/procA1.1/plot/bootstrapping_ticdim/bt_80_" + str(j) +"_msm.pkl",'rb'))
    feature = []
    for file in files:
        feature.extend(np.load(file)[:,[26,89,31,69]])

    feature = np.vstack(feature)
    boolean_list1 = feature[:,0] < threshold
    feat1 = list(map(int, boolean_list1))
    boolean_list2 = feature[:,1] < threshold
    feat2 = list(map(int, boolean_list2))
    boolean_list3 = feature[:,2] < threshold
    feat3 = list(map(int, boolean_list3))
    boolean_list4 = feature[:,3] < threshold
    feat4 = list(map(int, boolean_list4))
     
    weights = np.concatenate(msm.trajectory_weights())
    feat1_weighted =  np.dot(feat1, weights) #distance between Cys3 and Thr7
    feat2_weighted =  np.dot(feat2, weights) #distance between Thr12 and Cys16
    feat3_weighted =  np.dot(feat3, weights) #distance between Cys3 and Thr12
    feat4_weighted =  np.dot(feat4, weights) #distance between Thr7 and Cys16

    df.loc[j] = [feat1_weighted, feat2_weighted, feat3_weighted, feat4_weighted]
"""

#pickle.dump(df, open('procA1.1_four_feat_bt80_200.pkl','wb'))
df = pickle.load(open('procA1.1_four_feat_bt80_200.pkl','rb'))

df_np = df.to_numpy()


fig, axs = plt.subplots(1,1,figsize=(10,7))
axs.boxplot(df_np)
plt.xticks([1, 2, 3, 4], ['C3-T7', 'T12-C16', 'C3-T12', 'T7-C16'])
plt.ylabel('Probability of ring formation',fontsize=24)
plt.xticks(fontsize=24)
plt.yticks([0, 0.1, 0.2, 0.3, 0.4,0.5,0.6],fontsize=20)
plt.ylim(0, 0.6)
plt.savefig('procA1.1_box_plot.png')    

"""
feat1_mean = np.mean(df['feat1'])
feat2_mean = np.mean(df['feat2'])
feat3_mean = np.mean(df['feat3'])
feat4_mean = np.mean(df['feat4'])
feat1_std = np.std(df['feat1'])
feat2_std = np.std(df['feat2'])
feat3_std = np.std(df['feat3'])
feat4_std = np.std(df['feat4'])
feat = ['Cys3-Thr7','Thr12-Cys16','Cys3-Thr12','Thr7-Cys16']
feat_means = [feat1_mean, feat2_mean, feat3_mean, feat4_mean]
feat_error = [feat1_std, feat2_std, feat3_std, feat4_std]
x_pos = np.arange(len(feat))

fig, ax = plt.subplots()
ax.bar(x_pos, feat_means, yerr=feat_error, align='center', alpha=0.5, ecolor='black', capsize=10)
ax.set_ylabel('Probability of ring formation')
ax.set_xticks(x_pos)
ax.set_xticklabels(feat)
plt.tight_layout()
plt.savefig('procA1.1_bar_plot.png')
"""




