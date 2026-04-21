### Plot density for entropy wave ###

import numpy as np
import sys, glob, psutil, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import multiprocessing as mp
import warnings

warnings.filterwarnings("ignore")

mpl.rcParams['figure.dpi'] = 120
mpl.rcParams['savefig.dpi'] = 120
mpl.rcParams['figure.autolayout'] = True
mpl.rcParams['figure.figsize'] = (6,6)
mpl.rcParams['axes.titlesize'] = 18
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['text.usetex'] = False
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams["font.serif"] = 'cmr10',
mpl.rcParams["font.monospace"] = 'Computer Modern Typewriter'
mpl.rcParams["mathtext.fontset"]= 'cm'
mpl.rcParams['axes.unicode_minus'] = False

# paths
dumpsdir = sys.argv[1]
plotsdir = sys.argv[2]
if not os.path.exists(plotsdir):
	os.makedirs(plotsdir)

# function to parallelize plotting
def run_parallel(function, dlist,	nthreads):
	pool = mp.Pool(nthreads)
	pool.map_async(function, dlist).get(720000)
	pool.close()
	pool.join()
 

def plot(dumpno):
    print("Plotting dump {:04d}".format(dumpno))

    # header info
    header = open(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)), 'r')
    firstline = header.readline()
    header.close()
    firstline = firstline.split()

    n1   = int(firstline[6])
    n2   = int(firstline[7])
    t    = float(firstline[18])

    t = '{:.3f}'.format(t)

    # load density and four-current
    prims = np.loadtxt(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)), skiprows=1)
    rho  = prims[:,0].reshape((n1,n2))
    
    # print(rho[:,0])

	# read grid file
    grid = np.loadtxt(os.path.join(dumpsdir,'grid'))
    x1 = grid[:,4].reshape((n1,n2))
    x2 = grid[:,5].reshape((n1,n2))
    
	# plot	
    fig = plt.figure()
    nrows = 2
    ncols = 1
    heights = [1,10]
    gs = gridspec.GridSpec(nrows=nrows, ncols=ncols, height_ratios=heights, figure=fig)
    
    # plotting parameters
    # Expected entropy-wave density range: rho in [1-0.01, 1+0.01].
    vmin_rho = np.log10(1.0 - 1.e-4)
    vmax_rho = np.log10(1.0 + 1.e-4)
    cmap_rho = 'turbo'
    shading = 'gouraud'

    ax0 = fig.add_subplot(gs[0,:])
    ax0.annotate('t= '+str(t), xy=(0.5,0.5), xycoords='axes fraction', va='center', ha='center', fontsize = 'x-large')
    ax0.axis('off')

    ax1 = fig.add_subplot(gs[1,0])
    rho_plot = ax1.pcolormesh(x1, x2, np.log10(rho), cmap=cmap_rho, vmin=vmin_rho, vmax=vmax_rho, shading=shading)
    ax1.set_xlabel('$x$')
    ax1.set_ylabel('$y$')
    ax1.set_title('Log$_{10}(\\rho$)')
    ax1.set_aspect('equal')
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(rho_plot, cax=cax)

    plt.savefig(os.path.join(plotsdir,'sound_wave_{0:04d}.png'.format(dumpno)))
    plt.close()

if __name__=="__main__":
    dstart = int(sorted(glob.glob(os.path.join(dumpsdir,'dump*')))[0][-4:])
    dend = int(sorted(glob.glob(os.path.join(dumpsdir,'dump*')))[-1][-4:])
    dlist = range(dstart,dend+1)

    ncores = psutil.cpu_count(logical=True)
    pad = 0.5
    nthreads = int(ncores*pad)
    run_parallel(plot,dlist,nthreads)
