### Plot convergence for MHD linear modes in 2D ###

import numpy as np
import sys, glob, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
import click

@click.command()
@click.option('-r', '--res', default="64,128,256,512", help="Resolutions", show_default=True)
def plot_convergence(res):
    RES = []
    RES = [int(x) for x in res.split(",")]
    
    NVARS = 1
    VARS  = ['rho']
    
    l1 = np.zeros([len(RES), NVARS])
    powerfits = np.zeros(NVARS)
    
    for r in range(len(RES)):
        # load initial dump as reference solution
        dfile_init = sorted(glob.glob(os.path.join('{}'.format(RES[r]), 'dumps', 'dump_*')))[0]
        prims_init = np.loadtxt(dfile_init, skiprows=1)

        # load grid and final dump
        gfile = os.path.join('{}'.format(RES[r]), 'dumps','grid')
        dfile = sorted(glob.glob(os.path.join('{}'.format(RES[r]), 'dumps', 'dump_*')))[-1]
        prims = np.loadtxt(dfile, skiprows=1)
        
        # compute error
        for p in range(NVARS):
            l1[r,p] = np.mean(np.fabs(prims[:,p].reshape((RES[r],RES[r])) \
                - prims_init[:,p].reshape((RES[r],RES[r]))))
            
    # MEASURE CONVERGENCE
    for p in range(NVARS):
        powerfits[p] = np.polyfit(np.log(RES), np.log(l1[:,p]), 1)[0]
              
    # plotting parameters
    mpl.rcParams['figure.dpi']  = 120
    mpl.rcParams['savefig.dpi'] = 120
    mpl.rcParams['figure.autolayout'] = True
    mpl.rcParams['figure.figsize'] = (7,7)
    mpl.rcParams['axes.titlesize'] = 18
    mpl.rcParams['axes.labelsize'] = 16
    mpl.rcParams['xtick.labelsize'] = 14
    mpl.rcParams['ytick.labelsize'] = 14
    mpl.rcParams['legend.fontsize'] = 14
    colors = ['indigo', 'goldenrod', 'darkgreen', 'crimson', 'xkcd:blue', 'xkcd:magenta', 'green', 'xkcd:yellowgreen']
            
    # plot convergence
    fig = plt.gcf()
    ax  = fig.gca()
    
    ax.set_title('iharm2d_v4 advection')
    for p in range(NVARS):
        ax.loglog(RES, l1[:,p], color=colors[p], marker='o', label=VARS[p])
            
    ax.loglog([RES[0], RES[-1]], 5*np.asarray([RES[0], RES[-1]])**(-2.), color='black', linestyle='dashed', label='$N^{-2}$')
    plt.xscale('log', base=2)
    ax.set_xlabel('Resolution')
    ax.set_ylabel('L1 norm')
    ax.legend()
    plt.savefig(os.path.join("advection_2d_convergence.png"))


if __name__=='__main__':
    plot_convergence()


