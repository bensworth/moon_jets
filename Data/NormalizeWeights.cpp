#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <new>
#include <stdlib.h>
#include "gmp.h"
#include "gmpxx.h"
using namespace std;

// To compile, must include  -L /opt/local/lib -I /opt/local/include -lgmpxx -lgmp to link in
// arbitrary precision libraries.  
int main(int argc, char* argv[])
{
	double minVal = 1e-3;
	if(argc > 1)
	{
		if(strcmp(argv[1],"-minval") == 0)
			{
				minVal = atof(argv[2]);
		}
	}

	// Load .csv of weights for all particle speeds and sizes into a vector of vectors.
	mp_bitcnt_t precision = 2048;  // Create variable to set high precision in each array element
	mpf_class A[768][393];
	ifstream getWeight("./Weights.csv");
	assert(getWeight.is_open());
	int i=0,j=0;
	while(!getWeight.eof())
	{
		// Define temporary variables, get next line in .csv array.
		string line, tempS;
		size_t pos;
		getline(getWeight,line);	
		j = 0;
		// Pick out numeric substrings between commas, convert to double and save in vector.
		while( (pos=line.find(",")) != string::npos )
		{	
			tempS = line.substr(0,pos);			// Get substring
			const char * c = tempS.c_str();		// Convert to char*
			A[i][j].set_prec(precision);		// Set array element precision
			A[i][j].set_str(c, 10);				// Save as mpf variable in array
			line.erase(0,pos+1);
			j +=1;
		}
		// Get last element of line, pass current line into vector of vectors. 
		const char * c = tempS.c_str();
		A[i][j].set_prec(precision);
		A[i][j].set_str(c, 10);
		i += 1;
	}
	getWeight.close();

	
	mpf_class tNorm[768];	  // Create array to store sum of velocities at each particle size
	int bottom;
	// We will normalize twice. The first time we normalize all nontrivial rows. The second, we
	// will delete normalized entries less than machine precision, 1e-15, and renormalize. 
	for(int i=0; i<768; i++)
	{
		tNorm[i].set_prec(precision); 		// Set array element precision.
		for(int j=392; j>=0; j--)
		{
			tNorm[i] = tNorm[i] + A[i][j];	// Sum over a row (all velocities).
		}
		if(tNorm[i].get_d() < 1e-50)		// Break if rows have become arbitrarily small.
		{
			bottom = i-1;
			break; 				
		}
		for(int j=0; j<393; j++)
		{
			A[i][j] = A[i][j] / tNorm[i]; 	// Normalize velocities. 
		}
	}
	for(int i=0; i<=bottom; i++)
	{
		tNorm[i] = 0; 						// Rezero norm array.
		for(int j=392; j>=0; j--)
		{
			if(A[i][j].get_d() < minVal) A[i][j]= 0;	// Treat < minval as zero.
			tNorm[i] = tNorm[i] + A[i][j];				// Sum over a row (all velocities).
		}
		for(int j=0; j<393; j++)
		{
			A[i][j] = A[i][j] / tNorm[i]; 	// Renormalize velocities. 
		}
	}

	// // cout << tNorm[0].get_str(EXP) << ", " << EXP << ", " << cout << tNorm[0].get_prec() << endl;
	stringstream filename;
	filename << "./Weights_Normalized2.csv";
	ofstream output(filename.str());
	for(int i=0; i<=bottom; i++)
	{
		for(int j=0; j<392; j++)
		{
			output << A[i][j] << ",";
		}
		output << A[i][392] << endl;
	}


	return 0;
}
