import numpy as np
import sys, glob, psutil, os
from scipy import optimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import multiprocessing as mp

# Configure plot params
mpl.rcParams['figure.dpi'] = 120
mpl.rcParams['savefig.dpi'] = 120
mpl.rcParams['figure.autolayout'] = True
mpl.rcParams['figure.figsize'] = (8,5)
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

# Global dictionaries to store (i) fluid dump (ii) grid (iii) analytic solution data
dump = {}
grid = {}
soln = {}

############### GEOMETRY FUNCTIONS ###############
# Compute gcov in BL from (r,th) read from grid file
def gcov_bl():
    grid['gcov_bl'] = np.zeros_like(grid['gcov'])

    DD = 1 - 2./grid['r'] + grid['a']**2/grid['r']**2
    mu = 1 + grid['a']**2 * np.cos(grid['th'])**2 / grid['r']**2

    grid['gcov_bl'][Ellipsis,0,0] = -(1 - 2./(grid['r'] * mu))
    grid['gcov_bl'][Ellipsis,0,3] = -2 * grid['a'] * np.sin(grid['th'])**2 / (grid['r'] * mu)
    grid['gcov_bl'][Ellipsis,3,0] = grid['gcov_bl'][Ellipsis,0,3]
    grid['gcov_bl'][Ellipsis,1,1] = mu / DD
    grid['gcov_bl'][Ellipsis,2,2] = grid['r']**2 * mu
    grid['gcov_bl'][Ellipsis,3,3] = grid['r']**2 * np.sin(grid['th'])**2 * (1 + grid['a']**2/grid['r']**2 \
                                    + 2 * grid['a']**2 * np.sin(grid['th'])**2 / (grid['r']**3 * mu))


# Compute gcov in KS from (r,th) read from grid file
def gcov_ks():
    grid['gcov_ks'] = np.zeros_like(grid['gcov'])
    sigma = grid['r']**2 + (grid['a']**2 * np.cos(grid['th'])**2)
    
    grid['gcov_ks'][Ellipsis,0,0] = -1 + 2*grid['r']/sigma
    grid['gcov_ks'][Ellipsis,0,1] = 2*grid['r']/sigma
    grid['gcov_ks'][Ellipsis,0,3] = -(2*grid['a']*grid['r']*np.sin(grid['th'])**2)/sigma
    grid['gcov_ks'][Ellipsis,1,0] = 2*grid['r']/sigma
    grid['gcov_ks'][Ellipsis,1,1] = 1 + 2*grid['r']/sigma
    grid['gcov_ks'][Ellipsis,1,3] = -grid['a']*np.sin(grid['th'])**2 * (1 + 2*grid['r']/sigma)
    grid['gcov_ks'][Ellipsis,2,2] = sigma
    grid['gcov_ks'][Ellipsis,3,0] = -(2*grid['a']*grid['r']*np.sin(grid['th'])**2)/sigma
    grid['gcov_ks'][Ellipsis,3,1] = -grid['a']*np.sin(grid['th'])**2 * (1 + 2*grid['r']/sigma)
    grid['gcov_ks'][Ellipsis,3,3] = np.sin(grid['th'])**2 * (sigma + grid['a']**2*np.sin(grid['th'])**2 * (1 + 2*grid['r']/sigma))


# Compute gcov in KS from gcon_ks
def gcon_ks():
    grid['gcon_ks'] = np.linalg.inv(grid['gcov_ks'])


# Compute transformation matrix from KS -> MKS (for covariant indices)
def dxdX_KS_to_MKS():
    dxdX = np.zeros((grid['n1'], grid['n2'], grid['ndim'], grid['ndim']), dtype=float)

    dxdX[Ellipsis,0,0] = dxdX[Ellipsis,3,3] = 1
    dxdX[Ellipsis,1,1] = np.exp(grid['x1'])
    dxdX[Ellipsis,2,2] = np.pi + (1 - grid['hslope']) * np.pi * np.cos(2 * np.pi * grid['x2'])
    
    return dxdX


# Compute transformation matrix from MKS -> KS (for covariant indices)
def dxdX_MKS_to_KS():
    return (np.linalg.inv(dxdX_KS_to_MKS()))


# Compute quantities manually from x^mu
def bl_coords_from_x(grid_temp):
    grid_temp['r']  = np.exp(grid_temp['x1'])
    grid_temp['th'] = np.pi * grid_temp['x2'] + ((1 - grid['hslope'])/2.) * np.sin(2*np.pi*grid_temp['x2'])


def gcov_ks_from_x(grid_temp):
    bl_coords_from_x(grid_temp)

    grid_temp['gcov_ks'] = np.zeros_like(grid['gcov'])
    sigma = grid_temp['r']**2 + (grid_temp['a']**2 * np.cos(grid_temp['th'])**2)
    
    grid_temp['gcov_ks'][Ellipsis,0,0] = -1 + 2*grid_temp['r']/sigma
    grid_temp['gcov_ks'][Ellipsis,0,1] = 2*grid_temp['r']/sigma
    grid_temp['gcov_ks'][Ellipsis,0,3] = -(2*grid_temp['a']*grid_temp['r']*np.sin(grid_temp['th'])**2)/sigma
    grid_temp['gcov_ks'][Ellipsis,1,0] = 2*grid_temp['r']/sigma
    grid_temp['gcov_ks'][Ellipsis,1,1] = 1 + 2*grid_temp['r']/sigma
    grid_temp['gcov_ks'][Ellipsis,1,3] = -grid_temp['a']*np.sin(grid_temp['th'])**2 * (1 + 2*grid_temp['r']/sigma)
    grid_temp['gcov_ks'][Ellipsis,2,2] = sigma
    grid_temp['gcov_ks'][Ellipsis,3,0] = -(2*grid_temp['a']*grid_temp['r']*np.sin(grid_temp['th'])**2)/sigma
    grid_temp['gcov_ks'][Ellipsis,3,1] = -grid_temp['a']*np.sin(grid_temp['th'])**2 * (1 + 2*grid_temp['r']/sigma)
    grid_temp['gcov_ks'][Ellipsis,3,3] = np.sin(grid_temp['th'])**2 * (sigma + grid_temp['a']**2*np.sin(grid_temp['th'])**2 * (1 + 2*grid_temp['r']/sigma))


def dxdX_KS_to_MKS_from_x(grid_temp):
    dxdX = np.zeros((grid['n1'], grid['n2'], grid['ndim'], grid['ndim']), dtype=float)

    dxdX[Ellipsis,0,0] = dxdX[Ellipsis,3,3] = 1
    dxdX[Ellipsis,1,1] = np.exp(grid_temp['x1'])
    dxdX[Ellipsis,2,2] = np.pi + (1 - grid['hslope']) * np.pi * np.cos(2 * np.pi * grid_temp['x2'])

    return dxdX


def dxdX_MKS_to_KS_from_x(grid_temp):
    dxdX = dxdX_KS_to_MKS_from_x(grid_temp)
    return np.linalg.inv(dxdX)


def gcov_from_x(grid_temp):
    gcov_ks_from_x(grid_temp)
    dxdX = dxdX_KS_to_MKS_from_x(grid_temp)

    grid_temp['gcov'] = np.einsum('ijbn,ijmb->ijmn', dxdX, \
                        np.einsum('ijam,ijab->ijmb', dxdX, grid_temp['gcov_ks']))

    grid_temp['gcon'] = np.linalg.inv(grid_temp['gcov'])


# Compute the Christoffel symbols in MKS
def conn_func(sigma, alpha, beta):
    delta = 1.e-5
    conn = np.zeros((grid['n1'], grid['n2'], grid['ndim'], grid['ndim'], grid['ndim']), dtype=float)
    tmp  = np.zeros_like(conn)

    x = np.zeros((grid['n1'], grid['n2'], grid['ndim']), dtype=float)
    x[Ellipsis,1] = grid['x1']
    x[Ellipsis,2] = grid['x2']
    x[Ellipsis,3] = grid['x3']

    grid_h = {}; grid_h['a'] = grid['a']
    grid_l = {}; grid_l['a'] = grid['a']

    for mu in range(grid['ndim']):
        xh = np.copy(x)
        xl = np.copy(x)
        xh[Ellipsis,mu] += delta
        xl[Ellipsis,mu] -= delta

        grid_h['x1'] = xh[Ellipsis,1]
        grid_h['x2'] = xh[Ellipsis,2]
        grid_l['x1'] = xl[Ellipsis,1]
        grid_l['x2'] = xl[Ellipsis,2]

        gcov_from_x(grid_h)
        gcov_from_x(grid_l)

        for lam in range(grid['ndim']):
            for nu in range(grid['ndim']):
                conn[Ellipsis,lam,nu,mu] = (grid_h['gcov'][Ellipsis,lam,nu] - grid_l['gcov'][Ellipsis,lam,nu]) \
                                            / (xh[Ellipsis,mu] - xl[Ellipsis,mu])

    for lam in range(grid['ndim']):
        for nu in range(grid['ndim']):
            for mu in range(grid['ndim']):
                tmp[Ellipsis,lam,nu,mu] = 0.5 * (conn[Ellipsis,nu,lam,mu] + conn[Ellipsis,mu,lam,nu] \
                - conn[Ellipsis,mu,nu,lam])

    for lam in range(grid['ndim']):
        for nu in range(grid['ndim']):
            for mu in range(grid['ndim']):
                conn[Ellipsis,lam,nu,mu] = 0
                for kap in range(grid['ndim']):
                    conn[Ellipsis,lam,nu,mu] += grid['gcon'][Ellipsis,lam,kap] * tmp[Ellipsis,kap,nu,mu]

    return conn[Ellipsis,sigma,alpha,beta]


############### READ DATA ###############
# Read dump and/or grid file
def load_data(dumpsdir, dumpno):
    
    # header info
    header = open(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)), 'r')
    firstline = header.readline()
    header.close()
    firstline = firstline.split() 
    dump['mdot'] = float(firstline[0])   
    dump['rc']   = float(firstline[1])
    dump['gam']  = float(firstline[11])
    dump['rEH']  = float(firstline[21])
    
    grid['n1']   = float(firstline[7])
    grid['n2']   = float(firstline[8])
    grid['ndim'] = float(firstline[18])
    
    grid['dx1']  = float(firstline[16])
    grid['dx2']  = float(firstline[17])
    
    grid['rEH_ind'] = np.argmin(np.fabs(grid['r'][:,0] - dump['rEH']) > 0.)
    
    # read grid file
    grid = np.loadtxt(os.path.join(dumpsdir, 'grid'))
    grid['r']     = grid[:,2].reshape((grid['n1'],grid['n2']))
    grid['th']    = grid[:,3].reshape((grid['n1'],grid['n2']))
    grid['x1']    = grid[:,4].reshape((grid['n1'],grid['n2']))
    grid['x2']    = grid[:,5].reshape((grid['n1'],grid['n2']))
    grid['gdet']  = grid[:,6].reshape((grid['n1'],grid['n2']))
    grid['lapse'] = grid[:,7].reshape((grid['n1'],grid['n2']))
    grid['gcon']  = grid[:,8:24].reshape((grid['n1'], grid['n2'], grid['ndim'], grid['ndim']))
    grid['gcov']  = grid[:,24:].reshape((grid['n1'], grid['n2'], grid['ndim'], grid['ndim']))

    grid['x3'] = np.zeros((grid['n1'],grid['n2']))

    grid['hslope'] = float(firstline[23])
    grid['a'] = float(firstline[24])


############### COMPUTE ANALYTIC IDEAL BONDI SOLUTION ###############
# Nonlinear expression to solve for T
def T_func(T, r, C3, C4, N):
    return (1 + (1 + N/2)*T)**2 * (1 - 2./r + (C4**2/(r**4 * T**N))) - C3

# Obtain primitives for Bondi problem
def get_prim():
    N    = 2./ (dump['gam'] - 1)
    rc   = dump['rc']
    mdot = dump['mdot']
    vc   = np.sqrt(1. / (2 * rc))
    csc  = np.sqrt(vc**2 / (1 - 3*vc**2))
    Tc   = 2*N*csc**2 / ((N + 2)*(2 - N*csc**2))
    C4   = Tc**(N/2)*vc*rc**2
    C3   = (1 + (1 + N/2)*Tc)**2 * (1 - 2./rc + vc**2)

    # Root find T
    T = np.zeros_like(grid['r'][:,0])
    for index, r in enumerate(grid['r'][:,0]):
        T0       = Tc
        sol      = optimize.root(T_func, [T0], args=(r, C3, C4, N))
        T[index] = sol.x[0]
        if (sol.success!=True):
            print("Not converged at r = {:.2f}", r)

    # Compute remaining fluid variables
    soln['T'] = T
    soln['v'] = -C4 / (T**(N/2) * grid['r'][:,0]**2)
    soln['K'] = (4*np.pi*C4 / mdot) ** (2./N)

    soln['rho'] = soln['K']**(-N/2) * T**(N/2)
    soln['u']   = (N/2) * soln['K']**(-N/2) * T**(N/2 + 1)

    soln['mdot'] = mdot
    soln['N']    = N
    soln['rc']   = rc
    
    # We have u^r in BL. We need to convert this to ucon in MKS
    # First compute u^t in BL
    ucon_bl = np.zeros((grid['n1'], grid['n2'], grid['ndim']), dtype=float)
    AA = grid['gcov_bl'][Ellipsis,0,0]
    BB = 2. * grid['gcov_bl'][Ellipsis,0,1]*soln['v'][:,None]
    CC = 1. + grid['gcov_bl'][Ellipsis,1,1]*soln['v'][:,None]**2
    
    discr = BB*BB - 4.*AA*CC
    ucon_bl[Ellipsis,0] = (-BB - np.sqrt(discr)) / (2.*AA)
    ucon_bl[Ellipsis,1] = soln['v'][:,None]

    # Convert ucon(Bl) to ucon(KS)
    dxdX = np.zeros((grid['n1'], grid['n2'], grid['ndim'], grid['ndim']), dtype=float)
    dxdX[Ellipsis,0,0] = dxdX[Ellipsis,1,1] = dxdX[Ellipsis,2,2] = dxdX[Ellipsis,3,3] = 1.
    dxdX[Ellipsis,0,1] = 2*grid['r'] / (grid['r']**2 - 2.*grid['r'] + grid['a']**2)
    dxdX[Ellipsis,3,1] = grid['a']/(grid['r']**2 - 2.*grid['r'] + grid['a']**2)

    ucon_ks = np.zeros((grid['n1'], grid['n2'], grid['ndim']), dtype=float)
    for mu in range(grid['ndim']):
        for nu in range(grid['ndim']):
            ucon_ks[Ellipsis,mu] += dxdX[Ellipsis,mu,nu] * ucon_bl[Ellipsis,nu]

    # Convert ucon(KS) to ucon(MKS)
    ucon_mks = np.zeros((grid['n1'], grid['n2'], grid['ndim']), dtype=float)
    dxdX = dxdX_MKS_to_KS()
    for mu in range(grid['ndim']):
        for nu in range(grid['ndim']):
            ucon_mks[Ellipsis,mu] += dxdX[Ellipsis,mu,nu] * ucon_ks[Ellipsis,nu]

    ucov_mks = np.einsum('ijmn,ijn->ijm', grid['gcov'], ucon_mks)

    # Compute velocity primitives
    utilde = np.zeros((grid['n1'], grid['n2'], 3), dtype=float)

    alpha = 1./np.sqrt(-grid['gcon'][Ellipsis,0,0])
    beta  = np.zeros((grid['n1'], grid['n2'], 3), dtype=float)
    beta[Ellipsis,0] = alpha * alpha * grid['gcon'][Ellipsis,0,1]
    beta[Ellipsis,1] = alpha * alpha * grid['gcon'][Ellipsis,0,2]
    beta[Ellipsis,2] = alpha * alpha * grid['gcon'][Ellipsis,0,3]
    gamma = ucon_mks[Ellipsis,0] * alpha

    utilde[Ellipsis,0] = ucon_mks[Ellipsis,1] + beta[Ellipsis,0]*gamma/alpha
    utilde[Ellipsis,1] = ucon_mks[Ellipsis,2] + beta[Ellipsis,1]*gamma/alpha
    utilde[Ellipsis,2] = ucon_mks[Ellipsis,3] + beta[Ellipsis,2]*gamma/alpha
    
    soln['utilde_r'] = utilde[Ellipsis,0,0]
    
    
############### MAIN IS MAIN ###############
if __name__=='__main__':
    dumpsdir = sys.argv[1]   
    plotsdir = sys.argv[2]
    
    load_data(dumpsdir, 0, True)
    gcov_bl()
    gcov_ks()
    gcon_ks()
    get_prim()