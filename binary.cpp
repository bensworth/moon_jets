#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
using namespace std;

template<typename T> std::ostream& binary_write(std::ostream& stream, const T& value)
{
    return stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template<typename T> std::istream& binary_read(std::istream& stream, const T& value)
{
    return stream.read((char*) &value, sizeof(T));
}

int main(int argc, char *argv[])
{

// -------------------------------------------------------------------------------------- //
// --------------------------------- Read in Collisions --------------------------------- //
// -------------------------------------------------------------------------------------- //
#if 1

	ifstream C_stream("C-5_r438_a0.5.coll", ios::binary);
	int dataId, version, plasmaModel, jetRef, jetId, sizeId, numCores,
		runMinutes, total, collided, rows=0, cols=0, bodyID, MC, charging;
	float inclination, vgas, RC;
	vector<vector<float> > Collisions;
	while(!C_stream.eof())
	{
		int ID = -1;
		binary_read(C_stream,ID);
		switch(ID) { 
			case 0: 
				binary_read(C_stream,dataId);
				break;
			case 1: 
				binary_read(C_stream,version);
				break;
			case 2: 
				binary_read(C_stream, plasmaModel);
				break;
			case 3: 
				binary_read(C_stream, jetRef);
				break;
			case 4:
				binary_read(C_stream, jetId);
				break;
			case 5: 
				binary_read(C_stream, sizeId);
				break;
			case 6: 
				binary_read(C_stream, numCores);
				break;
			case 7:
				binary_read(C_stream, runMinutes);
				break;
			case 8: 
				binary_read(C_stream, total);
				break;
			case 9: 
				binary_read(C_stream, collided);		
				break;
			case 10: 
				binary_read(C_stream, inclination);
				break;
			case 11: 
				binary_read(C_stream, rows);
				binary_read(C_stream, cols);
				break;
			case 12: 
				if(rows==0 && cols==0) {
					cout << "Error -- Input size not defined. " << endl;
				}
				else {
					for(int i=0; i<rows; i++)
					{
						vector<float> temp(cols);
						for(int j=0; j<cols; j++)
						{
							float val; 
							binary_read(C_stream, val);
							temp[j] = val;
						}
						Collisions.push_back(temp);
					}
				}
				break;
			case 18:
				binary_read(C_stream, bodyID);
				break;
			case 19:
				binary_read(C_stream, MC);
				break;
			case 20:
				binary_read(C_stream, vgas);
				break;
			case 21:
				binary_read(C_stream, RC);
				break;
			case 23: 
				binary_read(C_stream, charging);
				break;
		}
	}
	C_stream.close();

	// Print results.
	for(int i=0; i<rows; i++)
	{
		for(int j=0; j<cols; j++)
		{
			cout << Collisions[i][j] << ", ";
		}
		cout <<endl;
	}	
	cout << "Data ID: " << dataId << "\nVersion: " << version << "\nPlasma Model: " << plasmaModel << "\nJet Ref. ID: " << jetRef << "\nJet ID: " << jetId << "\nSize ID: " << sizeId << "\nJet Inclination: " << inclination << "\nFraction collisions: " << collided << "/" << total << endl;
	cout << "Saved size: " << rows << " x " << cols << endl;
	cout << "Actual size: " << Collisions.size() << " x " << Collisions[0].size() << endl;

#endif

// ------------------------------------------------------------------------------------- //
// ---------------------------------- Read in Density ---------------------------------- //
// ------------------------------------------------------------------------------------- //
#if 0

	ifstream D_stream("D424_r381_a0.5.dens", ios::binary);
	// ifstream D_stream("D397_r381_a0.dens", ios::binary);
	int dataId, version, plasmaModel, jetRef, jetId, sizeId, numCores,
		runMinutes, size_x, size_y, size_z, min_x, min_y, min_z,
		angDistId, numAzimuth, N=0, bodyID, MC;
	float inclination, dx, critRad, gasSpeed;
	unordered_map<int,float> density; 
	vector<int> indices; 
	while(!D_stream.eof())
	{
		int ID = -1;
		binary_read(D_stream,ID);
		// cout << ID << endl;
		switch(ID) { 
			case 0: 
				binary_read(D_stream,dataId);
				// cout << "dataid " << dataId << endl;
				break;
			case 1: 
				binary_read(D_stream,version);
				cout << "version " << version << endl;
				break;
			case 2: 
				binary_read(D_stream, plasmaModel);
				break;
			case 3: 
				binary_read(D_stream, jetRef);
				cout << "jetRef " << jetRef << endl;
				break;
			case 4:
				binary_read(D_stream, jetId);
				cout << "id " << jetId << endl;
				break;
			case 5: 
				binary_read(D_stream, sizeId);
				cout << "sizeid " << sizeId << endl;
				break;
			case 6: 
				binary_read(D_stream, numCores);
				break;
			case 7:
				binary_read(D_stream, runMinutes);
				break;
			case 8: 
				binary_read(D_stream, size_x);
				binary_read(D_stream, size_y);
				binary_read(D_stream, size_z);
				break;
			case 22: 
				binary_read(D_stream, min_x);
				binary_read(D_stream, min_y);
				binary_read(D_stream, min_z);			
				break;
			case 10: 
				binary_read(D_stream, dx);
				break;
			case 11:
				binary_read(D_stream, angDistId);
				break;
			case 12:
				binary_read(D_stream, N);
				break;
			case 25: 
				if(N == 0) {
					cout << "Error -- N not defined. " << endl;
				}
				else {
					for(int i=0; i<N; i++) {
						long int temp; 
						binary_read(D_stream, temp);
						indices.push_back(temp);

						int z_ind = floor(temp / (size_x * size_y) ); 
						int y_ind = floor( (temp - z_ind * size_x * size_y) / size_x );
						int x_ind = temp - z_ind * size_x * size_y - y_ind * size_x; 

						// cout << "(" << x_ind << "," << y_ind << "," << z_ind << ") - " << endl;
					}
				}
				break;
			case 15: 
				if(N == 0) {
					cout << "Error -- N not defined. " << endl;
				}
				else {
					for(int i=0; i<N; i++) {
						float temp; 
						binary_read(D_stream, temp);
						density[indices[i]] += temp;
						// cout << density[indices[i]] << endl;
					}
				}
				break;
			case 16: 
				binary_read(D_stream, inclination);
				break;
			case 17: 
				binary_read(D_stream, numAzimuth);
				break;
			case 18: 
				binary_read(D_stream, bodyID);
				cout << "body " << bodyID << endl;
				break;
			case 19: 
				binary_read(D_stream, MC);
				break;
			case 20: 
				binary_read(D_stream, gasSpeed);
				break;
			case 21: 
				binary_read(D_stream, critRad);
				break;
		}
	}
	D_stream.close();
	
	cout << min_x << ", " << min_y << ", " << min_z << endl;
	cout << size_x << ", " << size_y << ", " << size_z << " - " << dx << endl;

	// // Print results. 
	// for(auto const & pair : density)
	// {	
	// 	int ind = pair.first; 
	// 	int z_ind = floor(ind / (size_x * size_y) ); 
	// 	int y_ind = floor( (ind - z_ind * size_x * size_y) / size_x );
	// 	int x_ind = ind - z_ind * size_x * size_y - y_ind * size_x; 

	// 	cout << "(" << x_ind << "," << y_ind << "," << z_ind << ") - " << pair.second << endl;
	// }	
	cout << N << " | " << density.size() << endl;

#endif


// ------------------------------------------------------------------------------------- //
// ----------------------------------- Write Density ----------------------------------- //
// ------------------------------------------------------------------------------------- // 
#if 0
	cent_x = 100; cent_y = 100; cent_z = 400; dx = 2.5;
	size_x = 200;	size_y = 200; size_z = 400;

	stringstream density_out;
	density_out << "./Jet39_fix/D" << jetId << "_r" << sizeId << "_a" << inclination << ".dens";
	ofstream dens_out(density_out.str(), ios::binary);

// Header Data. 
	int tmp;
	// Data ID
	binary_write(dens_out,tmp=0);		
	binary_write(dens_out,dataId);
	// Version
	binary_write(dens_out,tmp=1);		
	binary_write(dens_out,version);
	// Plasma model
	binary_write(dens_out,tmp=2);		
	binary_write(dens_out,plasmaModel);
	// Jet parameter reference
	binary_write(dens_out,tmp=3);		
	binary_write(dens_out,jetRef);
	// Jet ID
	binary_write(dens_out,tmp=4);		
	binary_write(dens_out,jetId); 
	// Size ID
	binary_write(dens_out,tmp=5);		
	binary_write(dens_out,sizeId);
	// Number of cores
	binary_write(dens_out,tmp=6);		
	binary_write(dens_out,numCores);
	// Grid data
	binary_write(dens_out,tmp=8);		
	binary_write(dens_out,size_x);
	binary_write(dens_out,size_y);
	binary_write(dens_out,size_z);
	// Center of grid
	binary_write(dens_out,tmp=9);		
	binary_write(dens_out,cent_x);
	binary_write(dens_out,cent_y);
	binary_write(dens_out,cent_z);
	// Dx
	binary_write(dens_out,tmp=10);		
	binary_write(dens_out,dx);
	// Length of data 
	binary_write(dens_out,tmp=12);		
	binary_write(dens_out,N);
	// Angle inclination
	binary_write(dens_out,tmp=16);		
	binary_write(dens_out,inclination);
	// Number of azimuth points for inclination
	binary_write(dens_out,tmp=17);		
	binary_write(dens_out,numAzimuth);
	// Density profile indices
	binary_write(dens_out,tmp=14);	
	for(auto const & pair : density)
	{
		int ind = pair.first;
		binary_write(dens_out,ind);	
	}
	// Density at respective indices
	binary_write(dens_out, tmp=15);	
	for(auto const & pair : density)
	{
		float val = pair.second;
		binary_write(dens_out,val);
	}	
	dens_out.close();	

#endif







	return 0;
}