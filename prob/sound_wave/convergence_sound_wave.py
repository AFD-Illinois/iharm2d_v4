### Plot convergence for sound wave test ###

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
    
    NVARS = 4
    VARS  = ['rho', 'u', 'u1', 'u2']
    
    k = 2.*np.pi
    amp = 1.e-4
    
    # Background state
    var_background = np.zeros(NVARS)
    var_background[0] = 1.
    var_background[1] = 1.
    var_background[2] = 0.
    var_background[3] = 0.
    
    # wavevector
    theta = np.pi/4.
    k *= np.sqrt(2.)
    k1 = k*np.cos(theta)
    k2 = k*np.sin(theta)
    
    # Relativistic sound speed
    # p0 = (gam-1)*u0, h0 = 1 + u0/rho0 + p0/rho0, cs^2 = gam*p0/(rho0*h0)
    gam_ad = 4./3.
    rho0   = var_background[0]
    u0     = var_background[1]
    p0     = (gam_ad - 1.) * u0
    h0     = 1. + u0/rho0 + p0/rho0
    cs     = np.sqrt(gam_ad * p0 / (rho0 * h0))   # = 2/sqrt(21) ~ 0.43643578
 
    # Phase velocity is simply cs * k_hat (a 3-velocity, no Lorentz factor)
    vp1 = cs * np.cos(theta)
    vp2 = cs * np.sin(theta)
 
    # print(f"cs = {cs:.12f}   vp1 = {vp1:.12f}   vp2 = {vp2:.12f}")
        
    # Perturbation
    dvar = np.zeros(NVARS)
    dvar[0] = 0.580429504019981   #  sqrt(63/187)
    dvar[1] = 0.773906005359975   #  sqrt(112/187)
    dvar[2] = 0.253320866973317   #  sqrt(12/187)
    dvar[2] *= np.cos(theta)
    dvar[3] = 0.253320866973317   #  sqrt(12/187)
    dvar[3] *= np.sin(theta)

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
            dvar_sol.append(dvar[p] * np.cos(k1*(x1 - vp1*t) + k2*(x2 - vp2*t)))
            l1[r,p] = np.mean(np.fabs(dvar_code[p] - dvar_sol[p]))
            
        # print(f"Res={RES[r]}, t={t}, x1 range=[{x1.min():.6f}, {x1.max():.6f}], x2 range=[{x2.min():.6f}, {x2.max():.6f}]")
        # print(f"  dvar_code: min={dvar_code[3].min():.6e}, max={dvar_code[3].max():.6e}")
        # print(f"  dvar_sol:  min={dvar_sol[3].min():.6e}, max={dvar_sol[3].max():.6e}")
        # print(f"  L1={l1[r,3]:.6e}")
                    
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
    
    ax.set_title('iharm2d_v4 sound_wave')
    for p in range(NVARS):
        ax.loglog(RES, l1[:,p], color=colors[p], marker='o', label=VARS[p])
            
    ax.loglog([RES[0], RES[-1]], 1e-2*np.asarray([RES[0], RES[-1]])**(-2.), color='black', linestyle='dashed', label='$N^{-2}$')
    plt.xscale('log', base=2)
    ax.set_xlabel('Resolution')
    ax.set_ylabel('L1 norm')
    ax.legend()
    plt.savefig(os.path.join("sound_wave_2d_convergence.png"))


if __name__=='__main__':
    plot_convergence()