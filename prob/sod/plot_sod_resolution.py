import glob
import multiprocessing as mp
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import psutil
import warnings

warnings.filterwarnings('ignore')

mpl.rcParams['figure.dpi'] = 120
mpl.rcParams['savefig.dpi'] = 120
mpl.rcParams['figure.autolayout'] = True
mpl.rcParams['figure.figsize'] = (10, 8)
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['text.usetex'] = False
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cmr10',
mpl.rcParams['font.monospace'] = 'Computer Modern Typewriter'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['axes.unicode_minus'] = False

RESOLUTIONS = [128, 256, 512, 4096]
COLORS = {128: 'goldenrod', 256: 'green', 512: 'crimson', 4096: 'black'}

script_dir = os.path.dirname(os.path.abspath(__file__))
dumps_dirs = {r: os.path.join(script_dir, str(r), 'dumps') for r in RESOLUTIONS}
plotsdir = os.path.join(script_dir, 'plots_resolution')
os.makedirs(plotsdir, exist_ok=True)


def dump_path(dumpsdir, dumpno):
    return os.path.join(dumpsdir, 'dump_0000{0:04d}'.format(dumpno))


def load_dump(dumpsdir, dumpno):
    path = dump_path(dumpsdir, dumpno)
    with open(path, 'r') as f:
        firstline = f.readline().split()

    tscale = float(firstline[0])
    n1 = int(firstline[6])
    n2 = int(firstline[7])
    gam = float(firstline[10])
    t = float(firstline[18]) * tscale

    prims = np.loadtxt(path, skiprows=1)
    rho = prims[:, 0].reshape((n1, n2))
    u = prims[:, 1].reshape((n1, n2))
    u1 = prims[:, 2].reshape((n1, n2))
    p = (gam - 1.0) * u

    grid = np.loadtxt(os.path.join(dumpsdir, 'grid'))
    x1 = grid[:, 4].reshape((n1, n2))

    x = x1[:, 0]
    rho_1d = rho[:, 0]
    p_1d = p[:, 0]
    u1_1d = u1[:, 0]
    u_1d = u[:, 0]

    isort = np.argsort(x)
    x = x[isort]
    rho_1d = rho_1d[isort]
    p_1d = p_1d[isort]
    u1_1d = u1_1d[isort]
    e_1d = u_1d[isort] / rho_1d[isort]

    p_1d /= tscale**2
    u1_1d /= tscale
    e_1d /= tscale**2

    return t, x, rho_1d, p_1d, u1_1d, e_1d


def plot(dumpno):
    print('Plotting dump {:04d}'.format(dumpno))

    fig, axs = plt.subplots(2, 2, sharex=True)
    t_label = None

    for res in RESOLUTIONS:
        ddir = dumps_dirs[res]
        path = dump_path(ddir, dumpno)
        if not os.path.exists(path):
            continue

        t, x, rho_1d, p_1d, u1_1d, e_1d = load_dump(ddir, dumpno)
        if t_label is None:
            t_label = t

        lw = 1.5 if res < 4096 else 1.0
        axs[0, 0].plot(x, rho_1d, lw=lw, color=COLORS[res], label=str(res))
        axs[0, 1].plot(x, p_1d, lw=lw, color=COLORS[res], label=str(res))
        axs[1, 0].plot(x, u1_1d, lw=lw, color=COLORS[res], label=str(res))
        axs[1, 1].plot(x, e_1d, lw=lw, color=COLORS[res], label=str(res))

    fig.suptitle('SOD 1D Profiles, t={:.3f}'.format(t_label if t_label is not None else 0), fontsize=16)

    axs[0, 0].set_ylabel(r'$\rho$')
    axs[0, 0].set_title('Density')
    axs[0, 1].set_ylabel(r'$p$')
    axs[0, 1].set_title('Pressure')
    axs[1, 0].set_xlabel('x1')
    axs[1, 0].set_ylabel(r'$v^1$')
    axs[1, 0].set_title('Velocity')
    axs[1, 1].set_xlabel('x1')
    axs[1, 1].set_ylabel(r'$u/\rho$')
    axs[1, 1].set_title('Specific internal energy')
    
    for ax in axs.flat:
        ax.set_xlim(0., 1.)
    axs[0, 0].set_ylim(0., 1.05)
    axs[0, 1].set_ylim(0., 0.42)
    axs[1, 0].set_ylim(-0.01, 0.61)
    axs[1, 1].set_ylim(0.7, 1.2)

    for ax in axs.flat:
        ax.grid(alpha=0.25)

    axs[0, 1].legend(title='N', fontsize=11, title_fontsize=11)

    plt.savefig(os.path.join(plotsdir, 'sod_{0:04d}.png'.format(dumpno)))
    plt.close(fig)


def run_parallel(function, dlist, nthreads):
    pool = mp.Pool(nthreads)
    pool.map_async(function, dlist).get(720000)
    pool.close()
    pool.join()


if __name__ == '__main__':
    # Use dump numbers present in the highest-resolution run as the reference set
    ref_dir = dumps_dirs[4096]
    dumpfiles = sorted(glob.glob(os.path.join(ref_dir, 'dump_*')))
    dlist = [int(os.path.basename(f).split('_')[-1]) for f in dumpfiles]

    ncores = psutil.cpu_count(logical=True)
    pad = 0.5
    nthreads = max(1, int(ncores * pad))
    run_parallel(plot, dlist, nthreads)
