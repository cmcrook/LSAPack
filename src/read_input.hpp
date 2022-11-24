#pragma once

#define NAME_LEN 256

#include <vector>

class read_input {

public:
	int dim;                     //dimension of simulation
	unsigned int seed;                    //Random seed
	int eventspercycle;          // # events per particle per cycle 
	int N;                        // # spheres
	double initialpf;               // initial packing fraction
	double maxpf;            // maximum packing fraction
	double temp;                      // initial temperature (temp=0 means v=0)
	double growthrate;
	double maxpressure;
	double maxcollisionrate;
	double maxSizeChange;
	double citer;
	std::vector<double> particle_sizes;         // ratio of sphere radii for bidisperse 
	std::vector<int> particle_counts;      // fraction of larger spheres
	std::vector<double> particle_masses;                 // ratio of sphere masses
	bool hardwalls;                   // false = periodic, true =hard walls
	char readfile[NAME_LEN];    // file with configuration; if new, creates new
	char distfile[NAME_LEN];    // file with configuration; if new, creates new
	char writefile[NAME_LEN];    // file to write configuration
	char datafile[NAME_LEN];       // file to write statistics

	int read(int argc, char* argv[]);

};
