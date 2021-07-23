#pragma once

#define NAME_LEN 256


#include <vector>

class read_input {

public:

  int eventspercycle;          // # events per particle per cycle 
  int N;                        // # spheres
  double initialpf;               // initial packing fraction
  double maxpf;            // maximum packing fraction
  double temp;                      // initial temperature (temp=0 means v=0)
  double growthrate;               
  double maxpressure;         
  double maxcollisionrate;    
  std::vector<double> particle_sizes;         // ratio of sphere radii for bidisperse 
  std::vector<double> particle_fractions;      // fraction of larger spheres
  std::vector<double> particle_masses;                 // ratio of sphere masses
  int hardwallBC;                   // =0 for periodic, =1 for hard walls
  char readfile[NAME_LEN];    // file with configuration; if new, creates new
  char writefile[NAME_LEN];    // file to write configuration
  char datafile[NAME_LEN];       // file to write statistics

  int read(int argc, char* argv[]);
 
};
