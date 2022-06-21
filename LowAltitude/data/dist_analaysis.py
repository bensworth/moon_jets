##############################################################
### 
###
##############################################################
import pdb
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

##############################################################
##############################################################
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

##############################################################
##############################################################
def getOneParticleFlux(hd_file, speed_size_fn, min_r=0, max_r=100):
	speeds = np.array(hd_file["speeds m per s"])
	radii = np.array(hd_file["radii m"])
	ns = hd_file.attrs["num speeds"]
	nr = hd_file.attrs["num radii"]
	if (nr != len(radii)):
		raise ValueError("Inconsistent HDF5 radii data!\n")
	if (ns != len(speeds)):
		raise ValueError("Inconsistent HDF5 speed data!\n")
	one_flux = np.zeros((hd_file.attrs["num altitudes"],hd_file.attrs["num inclinations"]))
	for i in range(0,ns-1):
		for j in range(0,nr-1):
			ds = speeds[i+1]-speeds[i]
			dr = radii[i+1]-radii[i]
			ss0 = speeds[i]
			ss1 = speeds[i+1]
			rr0 = radii[j]
			rr1 = radii[j+1]
			# Trapezoid quadrature rule, with heaviside function H(r-r0)
			if rr0 >= min_r and rr0 <= max_r:
				one_flux += speed_size_fn(rr0,ss0) * ds * dr * \
					np.array(hd_file["flux profiles"]["flux_r"+str(j)+"_s"+str(i)])
				one_flux += speed_size_fn(rr0,ss1) * ds * dr * \
					np.array(hd_file["flux profiles"]["flux_r"+str(j)+"_s"+str(i+1)])
			if rr1 >= min_r and rr1 <= max_r:
				one_flux += speed_size_fn(rr1,ss0) * ds * dr * \
					np.array(hd_file["flux profiles"]["flux_r"+str(j+1)+"_s"+str(i)])
				one_flux += speed_size_fn(rr1,ss1) * ds * dr * \
					np.array(hd_file["flux profiles"]["flux_r"+str(j+1)+"_s"+str(i+1)])
	one_flux *= 0.25
	return one_flux

def getOneParticleMassFlux(hd_file, speed_size_fn, rho, min_r=0, max_r=100):
	speeds = np.array(hd_file["speeds m per s"])
	radii = np.array(hd_file["radii m"])
	ns = hd_file.attrs["num speeds"]
	nr = hd_file.attrs["num radii"]
	if (nr != len(radii)):
		raise ValueError("Inconsistent HDF5 radii data!\n")
	if (ns != len(speeds)):
		raise ValueError("Inconsistent HDF5 speed data!\n")
	one_flux = np.zeros((hd_file.attrs["num altitudes"],hd_file.attrs["num inclinations"]))
	for i in range(0,ns-1):
		for j in range(0,nr-1):
			ds = speeds[i+1]-speeds[i]
			dr = radii[i+1]-radii[i]
			ss0 = speeds[i]
			ss1 = speeds[i+1]
			rr0 = radii[j]
			rr1 = radii[j+1]

			# Trapezoid quadrature rule, with heaviside function H(r-r0)
			if rr0 >= min_r and rr0 <= max_r:
				one_flux += speed_size_fn(rr0,ss0) * ds * dr * rr0 * rr0 * rr0 * \
					np.array(hd_file["flux profiles"]["flux_r"+str(j)+"_s"+str(i)])
				one_flux += speed_size_fn(rr0,ss1) * ds * dr * rr0 * rr0 * rr0 * \
					np.array(hd_file["flux profiles"]["flux_r"+str(j)+"_s"+str(i+1)])
			if rr1 >= min_r and rr1 <= max_r:
				one_flux += speed_size_fn(rr1,ss0) * ds * dr * rr1 * rr1 * rr1 * \
					np.array(hd_file["flux profiles"]["flux_r"+str(j+1)+"_s"+str(i)])
				one_flux += speed_size_fn(rr1,ss1) * ds * dr * rr1 * rr1 * rr1 * \
					np.array(hd_file["flux profiles"]["flux_r"+str(j+1)+"_s"+str(i+1)])
	# one_flux *= (2.0*rho/3.0) # BUG : this was initial implementation, think is wrong!
	one_flux *= (np.pi*rho/3.0)
	return one_flux

##############################################################
################### Distribution parameters ##################
##############################################################
# Power law size distribution and average volume of particle 
rmin = 0.5 * 1e-6 	# Minimum particle radius in m
rmax = 15 * 1e-6  	# Maximum particle radius in m
alpha = 3.1
Calpha =  (alpha-1) / (rmin**(1-alpha) - rmax**(1-alpha))
psize = lambda r: Calpha * r**(-alpha)
if alpha != 4.0:
	vav = ( 4*np.pi*(alpha-1.0)*(rmin**(4.0-alpha)-rmax**(4.0-alpha)) ) / \
		( 3.0*(alpha-4.0)*(rmin**(1.0-alpha)-rmax**(1.0-alpha)) )
else:
	vav = ( 4*np.pi*(alpha-1.0)*(np.log(rmax) - np.log(rmin)) ) / \
		( 3.0*(rmin**(1.0-alpha)-rmax**(1.0-alpha)) )

# Size-dependent speed distribution
rc = 0.8 * 1e-6		# Critical radius in um
vgas = 700			# gas velocity in m/s
pspeed = lambda v, r: (1.0 + r/rc)*(r/rc)*(v/(vgas*vgas))*(1-v/vgas)**(r/rc - 1.0)

# Combined speed-size distribution
pspeed_size = lambda r, v: pspeed(v,r) * psize(r)

##############################################################
################### Collection parameters ####################
##############################################################
# min_r = 0.5*1e-6		# Minimum particle radius to collect
min_r = 3.0*1e-6		# Maximum particle radius to collect
max_r = None
det_size = 25			# Detector size in cm
mass_prod = 25.0/80.0 	# Mass production kg/s = 25kg/s / 80 jets
max_alt = 1000			# Max altitude for cone heat map plots
# alt_set = [500, 1000, 2000, 3000]
alt_set = [100,200,300,400,500]  	# Altitudes to plot result at
show = False				# Show plot <-> show=True, save plot <-> show=False
angular_dist = "uni" 	# Options "uni" and "cos2"
mass_flux = True		# Plot results in mass flux. If false, plot particle flux


##############################################################
###################### Construct plots #######################
##############################################################
if angular_dist == "uni":
	hh = h5py.File('./EncFluxData_uni.hdf5', 'r')
else:
	hh = h5py.File('./EncFluxData.hdf5', 'r')
rho = 916.0 		  	# Density water ice in kg/m^3
sens_area = (det_size/100.0)**2  	# Sensitive area in m^2
scaling = sens_area * mass_prod / (rho*vav)
scaling *= 1000.0 		# Convert to g/s from kg/s
scaling *= 60.0 		# Convert to g/minute
if mass_flux:
	if max_r is not None:
		flux = getOneParticleMassFlux(hh, pspeed_size, rho, min_r, max_r)
		lbl = "Grams/minute, " + f'{min_r*1e6:.1f}' + " $<$ r $<$ " + f'{max_r*1e6:.1f}' + "$\mu$m"
		params = "rmin" + f'{min_r*1e6:.1f}' + "_rmax" + f'{max_r*1e6:.1f}' + \
			"_rc" + f'{rc*1e6:.1f}' + "_alpha" + str(alpha)
	else:
		flux = getOneParticleMassFlux(hh, pspeed_size, rho, min_r)
		lbl = "Grams/minute, r $>$ " + f'{min_r*1e6:.1f}' + "um"
		params = "rmin" + f'{min_r*1e6:.1f}' + \
			"_rc" + f'{rc*1e6:.1f}' + "_alpha" + str(alpha)
else:
	if max_r is not None:
		flux = getOneParticleFlux(hh, pspeed_size, min_r, max_r)
		lbl = "Particles/minute, " + f'{min_r*1e6:.1f}' + " $<$ r $<$ " + f'{max_r*1e6:.1f}' + "$\mu$m"
		params = "rmin" + f'{min_r*1e6:.1f}' + "_rmax" + f'{max_r*1e6:.1f}' + \
			"_rc" + f'{rc*1e6:.1f}' + "_alpha" + str(alpha)
	else:
		flux = getOneParticleFlux(hh, pspeed_size, min_r)
		lbl = "Particles/minute, r $>$ " + f'{min_r*1e6:.1f}' + "um"
		params = "rmin" + f'{min_r*1e6:.1f}' + \
			"_rc" + f'{rc*1e6:.1f}' + "_alpha" + str(alpha)
if angular_dist == "uni":
	params += "_uni.pdf"
else:
	params += ".pdf"
flux *= scaling
ttl = str(det_size)+"cm x "+str(det_size)+"cm detector, M$^+=$ " + str(mass_prod) + " kg/s"

##############################################################
nr = flux.shape[0]
nphi = flux.shape[1]
alt = np.array(hh['grid altitudes'])
alt_inds = np.where(alt < max_alt)[0]
alt = alt[alt_inds]
# Symmetrize flux for better viewing
incr = np.pi/180.0 * (90 - np.array(hh['grid inclinations']))
incl = np.pi/180.0 * (90 + np.array(hh['grid inclinations']))
aar, iir = np.meshgrid(alt, incr)
aal, iil = np.meshgrid(alt, incl)

# Plot flux cone for particle flux
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(projection='polar')
ax.set_thetamin(75)
ax.set_thetamax(105)
# Heat map
if mass_flux:
	pc = ax.pcolormesh(iir, aar, np.log(flux[alt_inds,:].T), shading='auto', vmin=-4, vmax=1)
	ax.pcolormesh(iil, aal, np.log(flux[alt_inds,:].T), shading='auto', vmin=-4, vmax=1)
	cbar = fig.colorbar(pc, ticks=[-4,-3,-2,-1,0,1], label=lbl)
	cbar.ax.set_yticklabels(['1e-4', '1e-3', '1e-2', '0.1', '1', '10'])  # vertically oriented colorbar)
else:
	pc = ax.pcolormesh(iir, aar, np.log(flux[alt_inds,:].T), shading='auto')
	ax.pcolormesh(iil, aal, np.log(flux[alt_inds,:].T), shading='auto')
# Contour lines
plt.contour(iil, aal, np.log10(flux[alt_inds,:].T), levels=[-4,-3,-2,-1,0,1], colors='k', linestyles='solid', linewidths=1)
plt.contour(iir, aar, np.log10(flux[alt_inds,:].T), levels=[-4,-3,-2,-1,0,1], colors='k', linestyles='solid', linewidths=1)
plt.title(ttl)
if show:
	plt.show()
else:
	if mass_flux:
		plt.savefig("./figures/mass_conePlot_"+params)
	else:
		plt.savefig("./figures/num_conePlot_"+params)

##############################################################
alt = np.array(hh['grid altitudes'])
na = len(alt_set)
inds = [(np.abs(alt - alt_set[i])).argmin() for i in range(0,na)]
sin_incs = np.sin((np.pi/180.0)*np.array(hh['grid inclinations']))

fig = plt.figure(figsize=(11,8))
ax = fig.add_subplot()
for i in range(0,na):
	dists = alt_set[i]*sin_incs
	color = next(ax._get_lines.prop_cycler)['color']
	plt.semilogy(dists, flux[inds[i],:], color=color)
	plt.semilogy(-dists, flux[inds[i],:], color=color, label=str(alt_set[i]) + "m hover")
if mass_flux:
	plt.ylim([1e-2,500])
	plt.xlim([-100,100])
plt.grid()
plt.title(ttl)
plt.ylabel(lbl)
plt.xlabel("Distance from plume center (m)")
plt.legend()
if show:
	plt.show()
else:
	if mass_flux:
		plt.savefig("./figures/mass_altPlot_"+params)
	else:
		plt.savefig("./figures/num_altPlot_"+params)

