import numpy as np
import h5py, sys, glob, psutil, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import multiprocessing as mp

# paths
dumpsdir = sys.argv[1]
outputdir = sys.argv[2]
if not os.path.exists(outputdir):
	os.makedirs(outputdir)

# function to parallelize plotting
def run_parallel(function, dlist,	nthreads):
	pool = mp.Pool(nthreads)
	pool.map_async(function, dlist).get(720000)
	pool.close()
	pool.join()

# function to generate 4-potential from magnetic field
def plotting_bfield_lines(ax, B1, B2, xp, zp, n1, n2, dx1, dx2, gdet, nlines=20):
	AJ_phi = np.zeros([n1,n2]) 
	for j in range(n2):
			for i in range(n1):
					AJ_phi[i,j] = (np.trapz(gdet[:i,j]*B2[:i,j],dx=dx1) - np.trapz(gdet[i,:j]*B1[i,:j],dx=dx2))
	AJ_phi -=AJ_phi.min()
	levels = np.linspace(0,AJ_phi.max(),nlines*2)
	ax.contour(xp, zp, AJ_phi, levels=levels, colors='k')

def plotting(dumpno):	
	# plotting parameters
	vmin = -5; vmax = 0
	domain = [0,50,-50,50]
	shading = 'gouraud'
	cmap = 'jet'
	bh = 'True'
	plt.clf()
	
	# header info
	header = open(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),'r')
	firstline = header.readline()
	header.close()
	firstline = firstline.split()
	
	madtype = int(firstline[0])
	rin = float(firstline[2])	
	rmax = float(firstline[3])	
	electrons = float(firstline[7])
	metric = firstline[9]
	n1 = int(firstline[11])
	n2 = int(firstline[12])
	
	# if electron heating was enabled
	if len(firstline) > 38:
		gam = float(firstline[20])
		dx1 = float(firstline[25])
		dx2 = float(firstline[26])
		ndim = int(firstline[27])
		if metric == 'FMKS':
			rEH = float(firstline[33])
			a = float(firstline[36])
			t = float(firstline[37])
		elif metric == 'MKS':
			rEH = float(firstline[30])
			a = float(firstline[33])
			t = float(firstline[34])

	# if electron heating was not enabled
	else:
		gam = float(firstline[15])
		dx1 = float(firstline[20])
		dx2 = float(firstline[21])
		ndim = int(firstline[22])	
		if metric == 'FMKS':
			rEH = float(firstline[28])
			a = float(firstline[31])
			t = float(firstline[32])
		elif metric == 'MKS':
			rEH = float(firstline[25])
			a = float(firstline[28])
			t = float(firstline[29])

	t = '{:.3f}'.format(t)
	
	# loading prims
	prims = np.loadtxt(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),skiprows=1)
	rho = prims[:,0].reshape((n1,n2))
	B1 = prims[:,5].reshape((n1,n2))
	B2 = prims[:,6].reshape((n1,n2))
	logrho = np.log10(rho)

	# reading grid file
	grid = np.loadtxt(os.path.join(dumpsdir,'grid'))
	x = grid[:,0].reshape((n1,n2))
	z = grid[:,1].reshape((n1,n2))
	r = grid[:,2].reshape((n1,n2))
	th = grid[:,3].reshape((n1,n2))
	x1 = grid[:,4].reshape((n1,n2))
	x2 = grid[:,5].reshape((n1,n2))
	gdet = grid[:,6].reshape((n1,n2))
	gcon = grid[:,8:24].reshape((n1, n2, ndim, ndim))
	gcov = grid[:,24:].reshape((n1, n2, ndim, ndim))
	
	xp = x 
	xp[:,0] = xp[:,-1] = 0
	zp = z
	fig = plt.figure(figsize=(5,10))
	heights = [1,10]
	gs = gridspec.GridSpec(nrows = 2, ncols = 1, height_ratios = heights, figure = fig)

	# plotting
	ax0 = fig.add_subplot(gs[0,:])
	ax0.annotate('t= '+str(t), xy=(0.5,0.5), xycoords='axes fraction', va='center', ha='center', fontsize = 'x-large')
	ax0.axis('off')
	
	ax1 = fig.add_subplot(gs[1,:])
	rho_plot = ax1.pcolormesh(xp, zp, logrho, cmap=cmap, vmin=vmin, vmax=vmax, shading=shading)
	plotting_bfield_lines(ax1, B1, B2, xp, zp, n1, n2, dx1, dx2, gdet, nlines = 5)
	ax1.set_xlabel('$x (GM/c^2)$')
	ax1.set_ylabel('$z (GM/c^2)$')
	ax1.set_xlim(domain[:2])
	ax1.set_ylim(domain[2:])
	ax1.set_title('Log($\\rho$)',fontsize='large')
	if bh:
					circle = plt.Circle((0,0), rEH, color='k')
					ax1.add_artist(circle)
	ax1.set_aspect('equal')
	divider = make_axes_locatable(ax1)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	plt.colorbar(rho_plot, cax=cax)
	
	plt.savefig(os.path.join(outputdir,'density_{0:04d}.png'.format(dumpno)))
	plt.close()

if __name__=="__main__":
	dstart = int(sorted(glob.glob(os.path.join(dumpsdir,'dump*')))[0][-4:])
	dend = int(sorted(glob.glob(os.path.join(dumpsdir,'dump*')))[-1][-4:])
	dlist = range(dstart,dend+1)

	ncores = psutil.cpu_count(logical=True)
	pad = 0.5
	nthreads = int(ncores*pad)
	run_parallel(plotting,dlist,nthreads)
