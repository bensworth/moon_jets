- My jet collision files are saved for a given angle, and correspond to the contribution of a single particle to impacts of Enceladus (kind of)
	+ For using Juergen’s speed weights, collisions = w / num_azimuth, where num_azimuth is the number of azimuthal angles simulated for this speed.
	+ For MonteCarlo, weights are just 1/num_azimuth, despite sampling num_speed difference speeds as well… So normalized not to one particle, but to num_speeds particles… This was only used for Europa —> Need to check on before running MC for Enceladus


- When jets are imported, they are weighted with an angular distribution, and then normalized by the surface area of a given patch for a 1 degree longitude x 1 degree latitude area to be in impacts per 1/m^2. This should mean that each surf file is the flux representative of a single particle. 

- When I plot collision rate, this file is loaded and multiplied by the production rate for a given particle size. It still corresponds to flux per m^2 (per second??)

- For the effective production rate of a jet, we scale by siting(jet)/sum(sitings) * num_jets. Not clear why we scale by num_jets / when it comes out.


Error with Jet 32 surface data


CHARGING
————————
	- Something seems to be wrong with size 348 for all Simon sims??






Enceladus Paper
————————————————
- Show we can reproduce within our constraints, CDA data on flybys with approximately
a single set of parameters on mass production and size distribution
- Confirms that speed distribution is relatively accurate
- Show which parts of fractures are active
- If we use bimodal size distribution, which largely agrees with having two power-law
like distributions for particles formed through nucleation, and larger particles, then
…
- Multiple size distributions is separate form Andy’s work, who used 1
- Always show error bars of CDA plots. Important to see that we are within 1 sigma
of measured data. 

	+ Look at spectral stuff from fitting code; Sascha thinks it should be commented
	out in there somewhere
- Handful of jets contributing to E21 flyby (1-10) are not pointed towards max impact 

Surface deposition
———————————————————
- Pattern shows that three body effects matter; would never get pattern from two-body
- Can we look at cratered terrain vs. not cratered terrain - shows where Enc. has been resurfaced
	+ Appears to be large smooth terrain between the wave structure of resurfacing
	+ Resurfaced at another time?