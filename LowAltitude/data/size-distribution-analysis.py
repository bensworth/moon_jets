import pdb
import os
import h5py
import numpy as np
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import fsolve, minimize_scalar, minimize
from scipy.signal import savgol_filter

fig_width_pt = 400.0
fontsize = 30
fontfamily = 'serif'
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*2      # height in inches
fig_size = [fig_width, fig_height]
params = {'backend': 'ps',
          'font.family': fontfamily,
          'font.serif':  'cm',
          'font.sans-serif': 'arial',
          'axes.labelsize': fontsize,
          'font.size': fontsize,
          'axes.titlesize': fontsize,
          'legend.fontsize': fontsize,
          'xtick.labelsize': fontsize,
          'ytick.labelsize': fontsize,
          'text.usetex': True,
          # 'figure.figsize': fig_size,
          'lines.linewidth': 2}
plt.rcParams.update(params)
lw = 4

########################################################
########################################################
########################################################

### Select flyby and scaling of binsize
CA7_time = 15
CA21_pretime = -8
CA21_posttime = 15
CDA_SCALE = 1.0
num_radii_per_bin = 20
rmax = 15.0 * CDA_SCALE * 1e-6
min_type = "all"    # Minimize over slope of distribution fit to each pair or HRD bins
# min_type = "opt"  # Minimize over alpha, optimized to fit power law at each altitude

# Initial guess
x0 = (4.5, 0.8*1e-6, 750)

# Bounds on alpha, rc, and vgas, respectively
bnds = ((3.0, 5.0), (2.5*1e-6, 5.0*1e-6), (600,1200))

########################################################
########################################################
########################################################

maxangle = 25
cdir = '/Users/southworth/Documents/Projects/Plumes/JetsCode/C-simulations/LowAltitude/data/'
hd = h5py.File(cdir + 'EncAltitudeData_max-angle' + str(maxangle) + '.hdf5', 'r')

# Functions to set speed and size distribution parameters to be used in called funcitons
alpha_fn = 0.0
rmin_fn = 0.0
rmax_fn = 0.0
def set_size_fn(a, r0, r1):
    global alpha_fn
    alpha_fn = a
    global rmin_fn
    rmin_fn = r0
    global rmax_fn
    rmax_fn = r1

rc_fn = 0.0
vgas_fn = 0.0
def set_speed_fn(rc, vgas):
    global rc_fn
    rc_fn = rc
    global vgas_fn
    vgas_fn = vgas

# Power law size distribution and average volume of particle 
def psize(r):
    return ((alpha_fn-1) / (rmin_fn**(1-alpha_fn) - rmax_fn**(1-alpha_fn))) * r**(-alpha_fn)

# Size-dependent speed distribution
def pspeed(r, v):
    return (1.0 + r/rc_fn)*(r/rc_fn)*(v/(vgas_fn*vgas_fn))*(1-v/vgas_fn)**(r/rc_fn - 1.0)

# Integrate speed-size distribution over speed w.r.t. speed_fn, that is \int d(r,s)p(s|r)
def integrateSpeed(distribution, radii, speeds, speed_fn):
    try: 
        nr = len(radii)
        quadrature = np.zeros((nr,))
    except:
        nr = 1
        quadrature = 0
    ns = len(speeds[speeds <= vgas_fn])
    for i in range(0,ns-1):
        ds = speeds[i+1]-speeds[i]
        ss0 = speeds[i]
        ss1 = speeds[i+1]
        quadrature += speed_fn(radii,ss0) * ds * distribution[:,i]
        quadrature += speed_fn(radii,ss1) * ds * distribution[:,i+1]

    quadrature *= 0.5
    return quadrature

# Integrate distribution of values in r over size_fn w.r.t. discrete indices as
# quadrature points. Returns scalar as approximate integral over the interval
# [ radii[indices[0]], radii[indices[1]] ]. dist_r should be 1d, already
# integrated over speed.
def integrateSize(dist_r, radii, indices, size_fn):
    quadrature = 0
    i0 = int(indices[0])
    i1 = int(indices[1])
    for i in range(i0,i1):
        dr = radii[i+1]-radii[i]
        rr0 = radii[i]
        rr1 = radii[i+1]
        # Trapezoid quadrature rule
#         print("quadrature: r=",radii[i],", dr=",dr,"p(ri)=",size_fn(rr0))  # DEBUG|
        quadrature += size_fn(rr0) * dr * dist_r[i]
        quadrature += size_fn(rr1) * dr * dist_r[i+1]

    quadrature *= 0.5
    return quadrature

### Load CDA data, solve for optimal alpha values
e21_sizes = CDA_SCALE*np.array([1.556,2.735,6.01]) * 1e-6
e7_sizes = CDA_SCALE*np.array([1.689,2.97,6.52]) * 1e-6

cdir = '/Users/southworth/Documents/Projects/Plumes/SizeDist'
E21_alpha_3 = np.loadtxt(cdir + '/cda_slopes/EN21_alpha_3.csv',delimiter=',');
E7_alpha_3 = np.loadtxt(cdir + '/cda_slopes/EN7_alpha_3.csv',delimiter=',');
e21_data = np.exp(E21_alpha_3[:,6:])
e7_data = np.exp(E7_alpha_3[:,6:])
e21_times = E21_alpha_3[:,0]
e7_times = E7_alpha_3[:,0]
e21_alts = E21_alpha_3[:,2] * 1e3
e7_alts = E7_alpha_3[:,2] * 1e3

# select only times < CA_time
# t_inds = np.where(e21_times < CA21_time)[0]
t_inds = np.where(np.logical_and(np.greater_equal(e21_times,CA21_pretime),\
    np.greater_equal(CA21_posttime,e21_times)))[0]
e21_alts = e21_alts[t_inds]
e21_times = e21_times[t_inds]
e21_data = e21_data[t_inds,:]

t_inds = np.where(np.abs(e7_times) < CA7_time)[0]
e7_alts = e7_alts[t_inds]
e7_times = e7_times[t_inds]
e7_data = e7_data[t_inds,:]

n21 = e21_data.shape[0]
n7 = e7_data.shape[0]

# Error functions to optimize alpha; for individual sensors, these are
# used to numerically perform a nonlinear solve, which should be exact.
# For all sensors (errE21 and errE7), we have 1 free parameter to fit
# two equations, so it is used as an optimization.
def M2err_e21(alpha, ind):
    cda = (e21_data[ind,1] + e21_data[ind,2]) / (e21_data[ind,0] + e21_data[ind,1] + e21_data[ind,2])
    model = (rmax**(1-alpha) - e21_sizes[1]**(1-alpha)) / \
        (rmax**(1-alpha) - e21_sizes[0]**(1-alpha))
    return (cda-model)*(cda-model) / (model*model)
def M3err_e21(alpha, ind):
    cda = e21_data[ind,2] / (e21_data[ind,0] + e21_data[ind,1] + e21_data[ind,2])
    model = (rmax**(1-alpha) - e21_sizes[2]**(1-alpha)) / \
        (rmax**(1-alpha) - e21_sizes[0]**(1-alpha))
    return (cda-model)*(cda-model) / (model*model)
def M2err_e7(alpha, ind):
    cda = (e7_data[ind,1] + e7_data[ind,2]) / (e7_data[ind,0] + e7_data[ind,1] + e7_data[ind,2])
    model = (rmax**(1-alpha) - e7_sizes[1]**(1-alpha)) / \
        (rmax**(1-alpha) - e7_sizes[0]**(1-alpha))
    return (cda-model)*(cda-model) / (model*model)
def M3err_e7(alpha, ind):
    cda = e7_data[ind,2] / (e7_data[ind,0] + e7_data[ind,1] + e7_data[ind,2])
    model = (rmax**(1-alpha) - e7_sizes[2]**(1-alpha)) / \
        (rmax**(1-alpha) - e7_sizes[0]**(1-alpha))
    return (cda-model)*(cda-model) / (model*model)
def errE21(alpha, ind):
    return M2err_e21(alpha, ind) + M3err_e21(alpha, ind)
def errE7(alpha, ind):
    return M2err_e7(alpha, ind) + M3err_e7(alpha, ind)

# Solve for alpha to reproduce each bin
alpha21 = np.zeros((n21,2))
alpha7 = np.zeros((n7,2))
for i in range(0,n21):
    alpha21[i,0] = fsolve(M2err_e21, x0=4, args=(i), xtol=1e-10)
    alpha21[i,1] = fsolve(M3err_e21, x0=4, args=(i), xtol=1e-10)
for i in range(0,n7):
    alpha7[i,0] = fsolve(M2err_e7, x0=4, args=(i), xtol=1e-10)
    alpha7[i,1] = fsolve(M3err_e7, x0=4, args=(i), xtol=1e-10)

# Minimize alpha to fit all three size bins
alpha21opt = np.zeros((n21,))
alpha7opt = np.zeros((n7,))
for i in range(0,n21):
    opt = minimize_scalar(errE21, args=(i), tol=1e-10, bounds=(2, 10), method='bounded')
    alpha21opt[i] = opt.x
for i in range(0,n7):
    opt = minimize_scalar(errE7, args=(i), tol=1e-10, bounds=(2, 10), method='bounded')
    alpha7opt[i] = opt.x

########################################################
########################################################
########################################################
# Data from HDF5 file
full_altitudes = np.array(hd['grid altitudes'])
full_densities = np.array(hd["altitude densities"])
speeds = np.array(hd['speeds m per s'])
num_speeds = len(speeds)

# Global parameters taht are changed within functions
all_params = {}
radii = []
densities = []
altitudes = []
cda_slopes = []
cda_optimal = []
bins = []
r0i = -1
r1i = -1
r2i = -1
r3i = -1

def set_flyby(flyby):
    global r0i, r1i, r2i, r3i
    global bins
    global radii
    global densities
    global altitudes
    global cda_slopes
    global cda_optimal

    if flyby == "e21":
        rmin = e21_sizes[0]
        radii1 = np.linspace(e21_sizes[0],e21_sizes[1],num_radii_per_bin)
        radii2 = np.linspace(e21_sizes[1],e21_sizes[2],num_radii_per_bin)
        radii3 = np.linspace(e21_sizes[2],rmax,num_radii_per_bin)
    elif flyby == "e7":
        rmin = e7_sizes[0]
        radii1 = np.linspace(e7_sizes[0],e7_sizes[1],num_radii_per_bin)
        radii2 = np.linspace(e7_sizes[1],e7_sizes[2],num_radii_per_bin)
        radii3 = np.linspace(e7_sizes[2],rmax,num_radii_per_bin)

    radii = np.concatenate((radii1[0:-1],radii2[0:-1],radii3))
    num_radii = len(radii)
    r0i = 0
    r1i = num_radii_per_bin-1
    r2i = 2*num_radii_per_bin-2
    r3i = 3*num_radii_per_bin-3
    # Get discrete bins aligned with simulation data close to E21/E7 bins,
    # ---> {[1.5,3),[3,6),[6,rmax)}
    bins = np.array([radii[r0i],radii[r1i],radii[r2i],radii[r3i]])

    # CDA altitudes from data (suppose we have this), get closest data
    # from simulated density profiles
    if flyby == "e7":
        data_altitudes = e7_alts
        cda_slopes = alpha7
        cda_optimal = alpha7opt
    else:
        data_altitudes = e21_alts
        cda_slopes = alpha21
        cda_optimal = alpha21opt

    na = len(data_altitudes)
    alt_inds = np.zeros((na,),dtype=int)
    for ai in range(0,na):
        alt_inds[ai] = np.argmin(np.abs(full_altitudes - data_altitudes[ai]))

    densities = full_densities[alt_inds][:]
    altitudes = full_altitudes[alt_inds][:]
        
def get_opt_slopes(alpha, rc, vgas):
    global r0i, r1i, r2i, r3i
    global bins
    global radii
    global densities
    global altitudes
    global cda_slopes
    global cda_optimal
    set_size_fn(alpha, radii[r0i], radii[r3i])
    set_speed_fn(rc, vgas)
    print(alpha_fn,",",rc_fn,",",vgas_fn)

    num_radii = len(radii)
    na = altitudes.shape[0]
    alpha_opt = np.zeros((na,))
    dist = np.zeros((na, num_radii))
    # Integrate distribution against speed distribution for each radius
    for rr in range(0,len(radii)):
        dist[:,rr] = integrateSpeed(densities, radii[rr], speeds, pspeed)

    # plt_interval = 10   # DEBUG : Plot PDF/CDF for every plt_interval altitudes
    for ai in range(0,na):

        # Integrate against size distribution in discrete bins in r
        pr1 = integrateSize(dist[ai,:], radii, [r0i,r3i], psize)
        pr2 = integrateSize(dist[ai,:], radii, [r1i,r3i], psize)
        pr3 = integrateSize(dist[ai,:], radii, [r2i,r3i], psize)
        pr2 /= pr1
        pr3 /= pr1
        # print(pr2, ",", pr3)

        # if ai%plt_interval == 0:
        #     plt.plot(radii*1e6,dist[ai,:]/np.linalg.norm(dist[ai,:]))

        # Error functions to optimize alpha; for individual sensors, these are
        # used to numerically perform a nonlinear solve, which should be exact.
        # For all sensors (errE21 and errE7), we have 1 free parameter to fit
        # three equations, so it is used as an optimization.
        def M2err(alpha):
            model = (bins[3]**(1-alpha) - bins[1]**(1-alpha)) / (bins[3]**(1-alpha) - bins[0]**(1-alpha))
            return (pr2-model)*(pr2-model) / (model*model)
        def M3err(alpha):
            model = (bins[3]**(1-alpha) - bins[2]**(1-alpha)) / (bins[3]**(1-alpha) - bins[0]**(1-alpha))
            return (pr3-model)*(pr3-model) / (model*model)
        def err(alpha):
            return M2err(alpha) + M3err(alpha)

        # Solve for alpha optimal alpha over all bins
        opt = minimize_scalar(err, options={'xatol':1e-10, 'disp':0}, bounds=(2, 10), method='bounded')
        alpha_opt[ai] = opt.x

    # plt.xlabel("radius (um)")
    # plt.show()
    # pdb.set_trace()

    return alpha_opt

def get_3slopes(alpha, rc, vgas):
    global r0i, r1i, r2i, r3i
    global bins
    global radii
    global densities
    global altitudes
    global cda_slopes
    global cda_optimal
    set_size_fn(alpha, radii[r0i], radii[r3i])
    set_speed_fn(rc, vgas)
    
    num_radii = len(radii)
    na = altitudes.shape[0]
    alpha_opt = np.zeros((na,2))
    dist = np.zeros((na, num_radii))
    # Integrate distribution against speed distribution for each radius
    for rr in range(0,len(radii)):
        dist[:,rr] = integrateSpeed(densities, radii[rr], speeds, pspeed)

    for ai in range(0,na):

        # Integrate against size distribution in discrete bins in r
        pr1 = integrateSize(dist[ai,:], radii, [r0i,r3i], psize)
        pr2 = integrateSize(dist[ai,:], radii, [r1i,r3i], psize)
        pr3 = integrateSize(dist[ai,:], radii, [r2i,r3i], psize)
        if pr1 == 0:
            pdb.set_trace()
        pr2 /= pr1
        pr3 /= pr1

        # Error functions to optimize alpha; for individual sensors, these are
        # used to numerically perform a nonlinear solve, which should be exact.
        # For all sensors (errE21 and errE7), we have 1 free parameter to fit
        # three equations, so it is used as an optimization.
        def M2err(alpha):
            model = (bins[3]**(1-alpha) - bins[1]**(1-alpha)) / (bins[3]**(1-alpha) - bins[0]**(1-alpha))
            return (pr2-model)*(pr2-model) / (model*model)
        def M3err(alpha):
            model = (bins[3]**(1-alpha) - bins[2]**(1-alpha)) / (bins[3]**(1-alpha) - bins[0]**(1-alpha))
            return (pr3-model)*(pr3-model) / (model*model)

        # Solve for alpha separate for each bin
        alpha_opt[ai,0] = fsolve(M2err, x0=4., xtol=1e-10)[0]
        alpha_opt[ai,1] = fsolve(M3err, x0=4., xtol=1e-10)[0]

    return alpha_opt

minerr = 1e10
params = []

# Function to minimize difference between least squares model
# predicted sloped and least squares CDA slopes.
def min_optimal(x):
    global params
    global minerr
    global radii
    global densities
    global altitudes
    global cda_slopes
    global cda_optimal
    model_slopes = get_opt_slopes(x[0],x[1],x[2])
    err = np.linalg.norm(model_slopes - cda_optimal)
    sc = np.linalg.norm(cda_optimal)
    print("rel err=",err/sc,", alpha=",x[0],", rc=", x[1]*1e6,", vgas=",x[2])
    if err < minerr:
        minerr = err
        params = x
    return err*err

# Function to minimize difference between three simulated slopes
# and three observed CDA slopes.
def min_slopes(x):
    global params
    global minerr
    global radii
    global densities
    global altitudes
    global cda_slopes
    global cda_optimal
    model_slopes = get_3slopes(x[0],x[1],x[2])
    err = np.linalg.norm(model_slopes - cda_slopes, ord='fro')
    sc = np.linalg.norm(cda_slopes)
    print("rel err=",err/sc,", alpha=",x[0],", rc=", x[1]*1e6,", vgas=",x[2])
    if err < minerr:
        minerr = err
        params = x
    return err*err


# Apply minimization to E21 and E7 separate
results_self_all = {}
results_self_opt = {}
for flyby in ["e21","e7"]:

    minerr = 1e10
    params = []
    set_flyby(flyby)

    # Solve minimization problem
    if min_type == "all":
        sol = minimize(min_slopes, x0=x0, options={'maxiter': 100}, bounds=bnds, method="Nelder-Mead")
    else:
        sol = minimize(min_optimal, x0=x0, options={'maxiter': 100}, bounds=bnds, method="Nelder-Mead")

    all_params[flyby] = params

    # Get slopes for optimal parameters
    results_self_all[flyby] = get_3slopes(all_params[flyby][0],all_params[flyby][1],all_params[flyby][2])
    results_self_opt[flyby] = get_opt_slopes(all_params[flyby][0],all_params[flyby][1],all_params[flyby][2])

# E7 based on E21 params
results_opp_all = {}
results_opp_opt = {}
set_flyby("e7")
results_opp_all["e7"] = get_3slopes(all_params["e21"][0],all_params["e21"][1],all_params["e21"][2])
results_opp_opt["e7"] = get_opt_slopes(all_params["e21"][0],all_params["e21"][1],all_params["e21"][2])

# E21 based on E7 params
set_flyby("e21")
results_opp_all["e21"] = get_3slopes(all_params["e7"][0],all_params["e7"][1],all_params["e7"][2])
results_opp_opt["e21"] = get_opt_slopes(all_params["e7"][0],all_params["e7"][1],all_params["e7"][2])

########################################################
########################################################
########################################################

# Plot using parameters optimized for flyby
cmap = plt.get_cmap("tab10")
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(20,7))
ax1.plot(e7_times, alpha7opt[:],color=cmap(0),label="HRD",linewidth=3)
ax1.fill_between(e7_times, np.min(alpha7[:,0:2],axis=1),np.max(alpha7[:,0:2],axis=1),color=cmap(0),alpha=0.35)

ax1.plot(e7_times, results_self_opt["e7"],color=cmap(1),label="Model",linewidth=3)
ax1.fill_between(e7_times, np.min(results_self_all["e7"][:,0:2],axis=1),np.max(results_self_all["e7"][:,0:2],axis=1),color=cmap(1),alpha=0.35)

ax1.set_ylim([3,5.5])
ax1.set_xlim([-CA7_time,CA7_time])
ax1.grid(True)
ax1.set_title("E7: {alpha="+"{:.2f}".format(all_params["e7"][0])+", rc="+"{:.2f}".format(all_params["e7"][1]*1e6)+"um, vgas="+"{:.0f}".format(all_params["e7"][2]))
ax1.set_ylabel("Binned size distribution slope")
ax1.set_xlabel("Time w.r.t. CA")
ax1.legend()

ax1r=ax1.twinx()
ax1r.plot(e7_times, e7_alts/1e3, 'k--')
ax1r.set_ylabel("Altitude (km) - -")

ax2.plot(e21_times, alpha21opt[:],color=cmap(0),label="HRD",linewidth=3)
ax2.fill_between(e21_times, np.min(alpha21[:,0:2],axis=1),np.max(alpha21[:,0:2],axis=1),color=cmap(0),alpha=0.35)

ax2.plot(e21_times, results_self_opt["e21"],color=cmap(1),label="Model",linewidth=3)
ax2.fill_between(e21_times, np.min(results_self_all["e21"][:,0:2],axis=1),np.max(results_self_all["e21"][:,0:2],axis=1),color=cmap(1),alpha=0.35)

ax2.set_ylim([3,5.5])
ax2.set_xlim([CA21_pretime ,CA21_posttime])
ax2.grid(True)
ax2.set_title("E21: {alpha="+"{:.2f}".format(all_params["e21"][0])+", rc="+"{:.2f}".format(all_params["e21"][1]*1e6)+"um, vgas="+"{:.0f}".format(all_params["e21"][2]))
ax2.set_ylabel("Binned size distribution slope")
ax2.set_xlabel("Time w.r.t. CA")
ax2.legend()

ax2r=ax2.twinx()
ax2r.plot(e21_times, e21_alts/1e3, 'k--')
ax2r.set_ylabel("Altitude (km) - -")

f.tight_layout()
# plt.savefig("size-distribution-optimal.png")
plt.savefig("size-distribution-optimal-hrd"+"{:.0f}".format(CDA_SCALE)+".png")
plt.clf()

# Plot using parameters optimized for other flyby
cmap = plt.get_cmap("tab10")
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(20,7))
ax1.plot(e7_times, alpha7opt[:],color=cmap(0),label="HRD",linewidth=3)
ax1.fill_between(e7_times, np.min(alpha7[:,0:2],axis=1),np.max(alpha7[:,0:2],axis=1),color=cmap(0),alpha=0.35)

ax1.plot(e7_times, results_opp_opt["e7"],color=cmap(1),label="Model",linewidth=3)
ax1.fill_between(e7_times, np.min(results_opp_all["e7"][:,0:2],axis=1),np.max(results_opp_all["e7"][:,0:2],axis=1),color=cmap(1),alpha=0.35)

ax1.set_ylim([3,5.5])
ax1.set_xlim([-CA7_time,CA7_time])
ax1.grid(True)
ax1.set_title("E7: {alpha="+"{:.2f}".format(all_params["e21"][0])+", rc="+"{:.2f}".format(all_params["e21"][1]*1e6)+"um, vgas="+"{:.0f}".format(all_params["e21"][2]))
ax1.set_ylabel("Binned size distribution slope")
ax1.set_xlabel("Time w.r.t. CA")
ax1.legend()

ax1r=ax1.twinx()
ax1r.plot(e7_times, e7_alts/1e3, 'k--')
ax1r.set_ylabel("Altitude (km) - -")

ax2.plot(e21_times, alpha21opt[:],color=cmap(0),label="HRD",linewidth=3)
ax2.fill_between(e21_times, np.min(alpha21[:,0:2],axis=1),np.max(alpha21[:,0:2],axis=1),color=cmap(0),alpha=0.35)

ax2.plot(e21_times, results_opp_opt["e21"],color=cmap(1),label="Model",linewidth=3)
ax2.fill_between(e21_times, np.min(results_opp_all["e21"][:,0:2],axis=1),np.max(results_opp_all["e21"][:,0:2],axis=1),color=cmap(1),alpha=0.35)

ax2.set_ylim([3,5.5])
ax2.set_xlim([CA21_pretime ,CA21_posttime])
ax2.grid(True)
ax2.set_title("E21: {alpha="+"{:.2f}".format(all_params["e7"][0])+", rc="+"{:.2f}".format(all_params["e7"][1]*1e6)+"um, vgas="+"{:.0f}".format(all_params["e7"][2]))
ax2.set_ylabel("Binned size distribution slope")
ax2.set_xlabel("Time w.r.t. CA")
ax2.legend()

ax2r=ax2.twinx()
ax2r.plot(e21_times, e21_alts/1e3, 'k--')
ax2r.set_ylabel("Altitude (km) - -")

f.tight_layout()
plt.savefig("size-distribution-swapped-hrd"+"{:.0f}".format(CDA_SCALE)+".png")
plt.clf()
pdb.set_trace()
