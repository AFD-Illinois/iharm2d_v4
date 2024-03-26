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
@click.option('-m', '--mode', default="fast", help="Eigenmode",\
    type=click.Choice(['entropy', 'slow', 'alfven', 'fast']), show_default=True)
def plot_convergence(res, mode):
    RES = []
    RES = [int(x) for x in res.split(",")]
    
    NVARS = 8
    VARS  = ['rho', 'u', 'u1', 'u2', 'u3', 'B1', 'B2', 'B3']
    
    amp = 1.e-4
    k1  = 2.*np.pi
    k2  = 2.*np.pi
    
    l1 = np.zeros([len(RES), NVARS])
    powerfits = np.zeros(NVARS)
    
    # Analytic solution
    # Background state
    var_background = np.zeros(NVARS)
    var_background[0] = 1.
    var_background[1] = 1.
    var_background[5] = 1.
    
    # pertubation
    dvar = np.zeros(NVARS)
    if mode == "entropy":
        dvar[0] = 1.
    elif mode == "slow":
        dvar[0] = 0.558104461559
        dvar[1] = 0.744139282078
        dvar[2] = -0.277124827421
        dvar[3] = 0.0630348927707
        dvar[5] = -0.164323721928
        dvar[6] = 0.164323721928
    elif mode == "alfven":
        dvar[4] = 0.480384461415
        dvar[7] = 0.877058019307
    elif mode == "fast":
        dvar[0] = 0.476395427447
        dvar[1] = 0.635193903263
        dvar[2] = -0.102965815319
        dvar[3] = -0.316873207561
        dvar[5] = 0.359559114174
        dvar[6] = -0.359559114174
    dvar *= amp
    
    for r in range(len(RES)):
        gfile = os.path.join('{}'.format(RES[r]), 'dumps','grid')
        dfile = sorted(glob.glob(os.path.join('{}'.format(RES[r]), 'dumps', 'dump_*')))[-1]
        
        grid = np.loadtxt(gfile)
        x1 = grid[:,0].reshape((RES[r],RES[r]))
        x2 = grid[:,1].reshape((RES[r],RES[r]))
        
        # loading prims and saving perturbations
        dvar_code = []
        prims = np.loadtxt(dfile, skiprows=1)
        for p in range(NVARS):
            dvar_code.append(prims[:,p].reshape((RES[r],RES[r])) - var_background[p])
            
        # eigenmode and error
        dvar_sol = []
        for p in range(NVARS):
            dvar_sol.append(np.real(dvar[p]*np.cos(k1*x1 + k2*x2)))
            l1[r,p] = np.mean(np.fabs(dvar_code[p] - dvar_sol[p]))
            
    # MEASURE CONVERGENCE
    for p in range(NVARS):
        if abs(dvar[p]) != 0.:
            powerfits[p] = np.polyfit(np.log(RES), np.log(l1[:,p]), 1)[0]
              
    # plotting parameters
    mpl.rcParams['figure.dpi']  = 120
    mpl.rcParams['savefig.dpi'] = 120
    mpl.rcParams['figure.autolayout'] = True
    mpl.rcParams['figure.figsize'] = (10,10)
    mpl.rcParams['axes.titlesize'] = 18
    mpl.rcParams['axes.labelsize'] = 16
    mpl.rcParams['xtick.labelsize'] = 14
    mpl.rcParams['ytick.labelsize'] = 14
    mpl.rcParams['legend.fontsize'] = 14
    colors = ['indigo', 'goldenrod', 'darkgreen', 'crimson', 'xkcd:blue', 'xkcd:magenta', 'green', 'xkcd:yellowgreen']
            
    # plot convergence
    fig = plt.gcf()
    ax  = fig.gca()
    
    ax.set_title('iharm2d_v4 MHD modes in 2D ({})'.format(mode))
    for p in range(NVARS):
        if abs(dvar[p]) != 0.:
            ax.loglog(RES, l1[:,p], color=colors[p], marker='o', label=VARS[p])
            
    ax.loglog([RES[0], RES[-1]], 10*amp*np.asarray([RES[0], RES[-1]])**(-2.), color='black', linestyle='dashed', label='$N^{-2}$')
    plt.xscale('log', base=2)
    ax.legend()
    plt.savefig(os.path.join("mhdmodes2d_{}_convergence.png".format(mode)))


if __name__=='__main__':
    plot_convergence()


