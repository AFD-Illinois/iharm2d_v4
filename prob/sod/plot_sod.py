import glob
import multiprocessing as mp
import os
import sys
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


dumpsdir = sys.argv[1]
plotsdir = sys.argv[2]
if not os.path.exists(plotsdir):
	os.makedirs(plotsdir)


def run_parallel(function, dlist, nthreads):
	pool = mp.Pool(nthreads)
	pool.map_async(function, dlist).get(720000)
	pool.close()
	pool.join()


def dump_path(dumpno):
	return os.path.join(dumpsdir, 'dump_0000{0:04d}'.format(dumpno))


def plot(dumpno):
	print('Plotting dump {:04d}'.format(dumpno))

	with open(dump_path(dumpno), 'r') as header:
		firstline = header.readline().split()

	tscale = float(firstline[0])
	n1 = int(firstline[6])
	n2 = int(firstline[7])
	gam = float(firstline[10])
	t = float(firstline[18]) * tscale

	prims = np.loadtxt(dump_path(dumpno), skiprows=1)
	rho = prims[:, 0].reshape((n1, n2))
	u = prims[:, 1].reshape((n1, n2))
	u1 = prims[:, 2].reshape((n1, n2))
	p = (gam - 1.0) * u

	grid = np.loadtxt(os.path.join(dumpsdir, 'grid'))
	x1 = grid[:, 4].reshape((n1, n2))

	# SOD is 1D along X1; use the only/first X2 index.
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
 
	# Rescale by tscale
	p_1d /= tscale**2
	u1_1d /= tscale
	e_1d /= tscale**2

	fig, axs = plt.subplots(2, 2, sharex=True)
	fig.suptitle('SOD 1D Profiles, t={:.3f}'.format(t), fontsize=16)

	axs[0, 0].plot(x, rho_1d, lw=1.5)
	axs[0, 0].set_ylabel('$\\rho$')
	axs[0, 0].set_title('Density')

	axs[0, 1].plot(x, p_1d, lw=1.5)
	axs[0, 1].set_ylabel('$p$')
	axs[0, 1].set_title('Pressure')

	axs[1, 0].plot(x, u1_1d, lw=1.5)
	axs[1, 0].set_xlabel('x1')
	axs[1, 0].set_ylabel('$v^1$')
	axs[1, 0].set_title('Velocity')

	axs[1, 1].plot(x, e_1d, lw=1.5)
	axs[1, 1].set_xlabel('x1')
	axs[1, 1].set_ylabel('$u/\\rho$')
	axs[1, 1].set_title('Specific internal energy')

	for ax in axs.flat:
		ax.grid(alpha=0.25)

	plt.savefig(os.path.join(plotsdir, 'sod_{0:04d}.png'.format(dumpno)))
	plt.close(fig)


if __name__ == '__main__':
	dumpfiles = sorted(glob.glob(os.path.join(dumpsdir, 'dump_*')))
	dlist = [int(os.path.basename(f).split('_')[-1]) for f in dumpfiles]

	ncores = psutil.cpu_count(logical=True)
	pad = 0.5
	nthreads = max(1, int(ncores * pad))
	run_parallel(plot, dlist, nthreads)
