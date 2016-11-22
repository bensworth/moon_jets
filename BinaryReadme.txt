Density Binary: 
	Id 0: Data ID			- 1 int
		+ 0 = Density for one inc. angle
		+ 1 = Aggregate density for particle size
	Id 1: Version 			— 1 int
		+ 1 = initial runs - probably not good data
		+ 2 = valid runs, saved jetID is 1 low
		+ 3 = everything appears to be right
		…
		+ 7 = Simon b-field,track particle charge
		      Error in density normazliation - factor of 100 too big
		+ 8 = current
	Id 2: Bfield model 		— 1 int
		+ 0 = No plasma model
		+ 1 = Connerney 
		+ 2 = Simon et al 
	Id 3: Jet reference		— 1 int
		+ 0 = Old, e.g. A1, B1,…
		+ 1 = Porco et al
	Id 4: Jet ID			- 1 int
	Id 5: Size ID			— 1 int 
	Id 6: Number of cores		— 1 int
	Id 7: Runtime minutes		— 1 int 
	Id 8: Dim			- 3 ints (size x, size y, size z)
	Id 9: Center (ONLY OLD DATA)	- 3 ints (size x, size y, size z) 
	Id 10: dx			- 1 float
	Id 11: Angular Distribution 	- 1 int
		+ 0 = uniform
		+ 1 = cos^2
	Id 12: Length of data		— 1 int
		+ Length: n
	Id 14: OLD - Density indices 	— n ints
	Id 25: NEW - Density indices 	— n long ints
	Id 15: Densities	 	— n floats
	Id 16: Inclination angle	- 1 float
	Id 17: Number azimuth		- 1 int
	Id 18: Body ID			- 1 int
		+ 1 = Enceladus
		+ 2 = Europa
	Id 19: Monte Carlo (T/F binary)	- 1 int
	Id 20: Gas speed (only for MC)  - 1 float 
	Id 21: Critical radius (MC)  	- 1 float 
	Id 22: Minimum grid values 	- 3 ints (min x, min y, min z) 
	Id 23: Particle Charging	- 1 int
		+ 1 = No charging
		+ 2 = Constant charge
		+ 3 = Full charging
	ID 24: Particle charge 		- n floats

Collision Binary: 
	Id 0: Data ID			- 1 int
	Id 1: Version 			— 1 int
	Id 2: Bfield model 		— 1 int 
	Id 3: Jet reference		— 1 int
	Id 4: Jet ID			- 1 int
	Id 5: Size ID			— 1 int 
	Id 6: Number of cores		— 1 int
	Id 7: Runtime minutes		— 1 int 
	Id 8: Launched particles	— 1 long int
	Id 9: Collided particles 	— 1 long int
	Id 10: Inclination angle	- 1 float
	Id 10: Inclination angle	- 1 float
	Id 11: Dimensions of data 	— 2 ints: rows, columns
	Id 12: Collision data		— (rows*columns) floats
	Id 18: Body ID			- 1 int
		+ 1 = Enceladus
		+ 2 = Europa
	Id 19: Monte Carlo (T/F binary)	- 1 int
	Id 20: Gas speed (only for MC)  - 1 float 
	Id 21: Critical radius (MC)  	- 1 float
	Id 23: Particle Charging	- 1 int
		+ 1 = No charging
		+ 2 = Constant charge
		+ 3 = Full charging

Density index to coordinates:

	int NUM_x = 200, NUM_y = 200, NUM_z = 600, DX = 10;
	z_ind = floor(ind / (NUM_x * NUM_y) );	
	y_ind = floor( (ind - z_ind * NUM_x * NUM_y) / NUM_x );
	x_ind = ind - z_ind * NUM_x * NUM_y - y_ind * NUM_x; 
	x_coor = DX*x_ind + .5*DX * (1. - NUM_x);
	y_coor = DX*y_ind + .5*DX * (1. - NUM_y);
	z_coor = DX*z_ind + .5*DX * (1. - NUM_z);

	- Note, the coordinates are Euclidean coordinates in the center of the grid.