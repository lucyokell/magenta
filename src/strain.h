//
//  MAGENTA
//  strain.h
//
//  Created: OJ on 13/01/2017
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Class describing a parasite strain.
//
// ---------------------------------------------------------------------------


#ifndef STRAIN_H
#define STRAIN_H

#include <vector>
#include <queue>
#include <bitset>
//using namespace std;

const int barcode_length = 24;
using barcode_t = std::bitset<barcode_length>;

class Strain {


private:

	static int s_strain_ID_generator;	// static member variable for initialising unique identifiers
	const int m_strain_ID;				// member variable ID	
	barcode_t* m_barcode_pointer;		// barcode sequence pointer  

	// Phenotype parameters to do here for resistance work
	// TODO:

public:

	// Default class constructor
	Strain();

	//


};

#endif
