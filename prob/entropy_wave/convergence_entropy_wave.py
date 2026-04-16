### Plot convergence for entropy wave test ###

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
    
    k = 2.*np.pi
    amp = 1.e-2
    
    # Background state
    var_background = np.zeros(NVARS)
    var_background[0] = 1.
    
    # wavevector
    theta = np.pi/4.
    k *= np.sqrt(2.)
    k1 = k*np.cos(theta)
    k2 = k*np.sin(theta)
    
    # Get wave velocity
    u10 = 0.1
    u1 = u10 * np.cos(theta)
    u2 = u10 * np.sin(theta)
    gamma = np.sqrt(1.0 + u1**2 + u2**2)  # √(1.01) ≈ 1.00499
    v1 = u1 / gamma
    v2 = u2 / gamma
    
    # Perturbation
    dvar = np.zeros(NVARS)
    dvar[0] = 1.
    
    dvar *= amp
    
    l1 = np.zeros([len(RES), NVARS])
    powerfits = np.zeros(NVARS)
    
    for r in range(len(RES)):
        gfile = os.path.join('{}'.format(RES[r]), 'dumps','grid')
        dfile = sorted(glob.glob(os.path.join('{}'.format(RES[r]), 'dumps', 'dump_*')))[-1]
        
        # get timestamp
        header = open(dfile, 'r')
        firstline = header.readline()
        header.close()
        firstline = firstline.split()
        t = float(firstline[18])
                                
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
            dvar_sol.append(dvar[p]*np.cos(k1*(x1 - v1*t) + k2*(x2 - v2*t)))
            l1[r,p] = np.mean(np.fabs(dvar_code[p] - dvar_sol[p]))
            
        # print(f"Res={RES[r]}, t={t}, x1 range=[{x1.min():.6f}, {x1.max():.6f}], x2 range=[{x2.min():.6f}, {x2.max():.6f}]")
        # print(f"  dvar_code: min={dvar_code[0].min():.6e}, max={dvar_code[0].max():.6e}")
        # print(f"  dvar_sol:  min={dvar_sol[0].min():.6e}, max={dvar_sol[0].max():.6e}")
        # print(f"  L1={l1[r,0]:.6e}")
                    
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
    
    ax.set_title('iharm2d_v4 entropy_wave')
    for p in range(NVARS):
        ax.loglog(RES, l1[:,p], color=colors[p], marker='o', label=VARS[p])
            
    ax.loglog([RES[0], RES[-1]], 1*np.asarray([RES[0], RES[-1]])**(-2.), color='black', linestyle='dashed', label='$N^{-2}$')
    plt.xscale('log', base=2)
    ax.set_xlabel('Resolution')
    ax.set_ylabel('L1 norm')
    ax.legend()
    plt.savefig(os.path.join("entropy_wave_2d_convergence.png"))


if __name__=='__main__':
    plot_convergence()