import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

import argparse
import glob
plt.rc('savefig', dpi=500)
matplotlib.rc('font',family='Helvetica-Normal',size=13)

name = 'procA3.3WT_cluster_150_ticdim_10'
data_name='procA3.3WT-four-dist-sub-helical'

def get_prob_density(data1, data2, bins=100, weights=None):

    if len(np.shape(data1)) > 1:
        data1 = np.asarray(data1)[:,0]
        data2 = np.asarray(data2)[:,0]

    hist, x_edges, y_edges = np.histogram2d(data1, data2, bins=100, weights=weights)

    prob_density = hist/np.sum(hist)
    x_coords = 0.5*(x_edges[:-1]+x_edges[1:]) 
    y_coords = 0.5*(y_edges[:-1]+y_edges[1:])

    return prob_density.T, x_coords, y_coords

def free_energy(data1, data2, T=300, weights=None, lims = None, max_energy=5, label1="Label Your Axes!", label2="Label Your Axes!", savename="fe.png", title = None):

    #compute free energy
    gas_constant = 0.00198588
    prob, x, y = get_prob_density(data1, data2, weights=weights)
    X, Y = np.meshgrid(x,y)
    free_energy = -gas_constant*T*np.log(prob)
    free_energy -= np.min(free_energy)
    plt.figure()
    fig, ax = plt.subplots()
    
    CS = plt.contourf(X, Y, free_energy, np.linspace(0, max_energy, max_energy*5+1), vmin=0.0, vmax=max_energy, cmap='jet')
    cbar = plt.colorbar(ticks=range(max_energy+1))
    cbar.set_label("Free Energy (kcal/mol)",size=16)
    cbar.ax.set_yticklabels(range(max_energy+1))
    cbar.ax.tick_params(labelsize=16)
    
    #plt.contour(X,Y,free_energy, linewidths=0.25, levels=[0,1,2,3,4,5], colors='black')
 
    plt.tick_params(axis='both',labelsize=12)
    plt.xlabel(label1)
    plt.ylabel(label2)
    #plt.ylim(1,9)
    #plt.xlim(0,1)
    plt.xlim(lims[0],lims[1])
    plt.ylim(lims[2],lims[3])
    plt.title(title)
    #axes = plt.axes()
    
    
    if lims[1] - lims[0] <= 1:
        ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
    elif lims[1] - lims[0] <= 2:
        ax.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    else:
        ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))

    if lims[3] - lims[2] <= 1:
        ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    elif lims[3] - lims[2] <= 2:
        ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    else:
        ax.yaxis.set_major_locator(plt.MultipleLocator(5.0))
    
    plt.grid(linestyle=":")
    #ax.spines['bottom'].set_linewidth(1)
    #ax.spines['left'].set_linewidth(1)
    #ax.spines['top'].set_linewidth(1)
    #ax.spines['right'].set_linewidth(1)
    #plt.gca().set_aspect(aspect=1.0, adjustable='box')
    fig.tight_layout()

    #plt.plot(4.9,2.7, 'ko', markersize=4, mec='w', mew=0.3)
    plt.plot(0.33,5.7, 'm*', markersize=16, mec='w', mew=0.3)


    plt.savefig(savename, transparent=False)


if __name__=="__main__":


    ax1_data = []
    ax2_data = []
    import pickle
    data = pickle.load(open(data_name +'.pkl', 'rb'), encoding='latin1')
    for i in range(len(data)):
        for j in range(len(data[i])):
             ax1_data.append(data[i][j][4])
             ax2_data.append(data[i][j][1]*10)
    
    print(np.shape(ax1_data))
    print(np.shape(ax2_data))

    
    weights_all = pickle.load(open('MSM-'+name+'-weights.pkl','rb'), encoding='latin1')
    weights = []
    for traj in weights_all:
        weights.extend(traj)
    weights = np.array(weights)
    print(np.shape(weights))
   
    free_energy(ax1_data, ax2_data,lims = (0,1,1,31), max_energy=5, label1='Helical content', label2= r'Distance between Thr11 and Cys21 ($\AA$)', savename=data_name+'-'+name+'-T11-C21.png', weights=weights, title = '')
