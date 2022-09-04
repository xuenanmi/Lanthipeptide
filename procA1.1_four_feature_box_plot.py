import numpy as np
import pandas as pd
import glob
import mdtraj as md
import matplotlib.pyplot as plt
import pickle


df = pd.DataFrame(columns = ['feat1','feat2','feat3','feat4'])
threshold = 0.45


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


pickle.dump(df, open('procA1.1_four_feat_bt80_200.pkl','wb'))
#df = pickle.load(open('procA1.1_four_feat_bt80_200.pkl','rb'))

df_np = df.to_numpy()


fig, axs = plt.subplots(1,1,figsize=(10,7))
axs.boxplot(df_np)
plt.xticks([1, 2, 3, 4], ['C3-T7', 'T12-C16', 'C3-T12', 'T7-C16'])
plt.ylabel('Probability of ring formation',fontsize=24)
plt.xticks(fontsize=24)
plt.yticks([0, 0.1, 0.2, 0.3, 0.4,0.5,0.6],fontsize=20)
plt.ylim(0, 0.6)
plt.savefig('procA1.1_box_plot.png')    






