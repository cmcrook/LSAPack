//===========================================================
//===========================================================
//===========================================================
//
//  Molecular dynamics simulation of hardspheres
//
//===========================================================
//===========================================================
//===========================================================

#define strncasecmp _strnicmp
#define strcasecmp _stricmp

#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <time.h>
#include <string.h>

#include "box.hpp"
#include "sphere.hpp"
#include "event.hpp"
#include "heap.hpp"
#include "read_input.hpp"


int main(int argc, char **argv)
{
  read_input input;
  int error = input.read(argc, argv);
  if (error) return error;

  double d, r;   // initial diameter and radius of spheres

  if(strcasecmp(input.readfile, "new")==0)
    input.readfile[0]=0;

  if (input.readfile[0]) // read in existing configuration
    {
      // read the header
      std::ifstream infile(input.readfile);
      if (!infile)
	{
	  std::cout << "error, can't open " << input.readfile  << std::endl;
	  exit(-1);
	}
      else
	{
	  int dim;
	  infile >> dim; infile.ignore(256, '\n');
	  if (dim != DIM)  // quit if dimensions don't match
	    {
	      std::cout << "error, dimensions don't match" << std::endl;
	      exit(-1);
	    }
	  infile.ignore(256, '\n');  // ignore the N 1 line
	  infile >> input.N; infile.ignore(256, '\n');
	  std::cout << "N = " << input.N << std::endl;
	  infile >> d; infile.ignore(256, '\n');
	  std::cout << "d = " << d << std::endl;
	  r = d/2.;
	  std::cout << "r = " << r << std::endl;
	}
    }
  else // create a new configuration
    {
      r = pow(input.initialpf*pow(SIZE, DIM)/(input.N*VOLUMESPHERE), 1.0/((double)(DIM)));
    }

  //Temporary radii distribution for testing
  //std::vector<double> particle_sizes = { 1.0 };
  //std::vector<double> particle_fractions = { 1.0 };
  //std::vector<double> particle_masses = {1.0 };

  std::vector<double> particle_sizes = { 0.02, 0.1, 0.2, 0.5, 1.0 };
  std::vector<double> particle_fractions = { 0.685, 0.2, 0.1, 0.01, 0.005 };
  std::vector<double> particle_masses(particle_sizes.size(), 1.0);

  box b(input.N, r, input.growthrate, input.maxpf, particle_sizes,
      particle_fractions, particle_masses, input.hardwallBC);
  
  std::cout << "ngrids = " << b.ngrids << std::endl;
  std::cout << "DIM = " << DIM << std::endl;

  if(input.readfile[0])
    {
      std::cout << "Reading in positions of spheres" << std::endl;
      b.RecreateSpheres(input.readfile, input.temp);
    }
  else 
    {
      std::cout << "Creating new positions of spheres" << std::endl;
      b.CreateSpheres(input.temp);
    } 
  
  std::ofstream output(input.datafile);
  output.precision(16);  

  std::cout << "Iteration " <<  "Collision Rate " << "Packing Fraction " << "Pressure " << std::endl;

  int iteration = 0;
  b.WriteConfiguration(input.writefile, iteration);
  while ((b.collisionrate < input.maxcollisionrate) && (b.pf < input.maxpf) && (b.pressure < input.maxpressure)) 
    {
      b.Process(input.eventspercycle*input.N);
      output << b.pf << " " << b.pressure << " " << 
	  b.collisionrate << " " << b.neventstot << " " << std::endl;

      b.Synchronize(true);
      std::cout << iteration << " " << b.collisionrate << " " << b.pf << "  " << b.pressure << std::endl;

      iteration++;
      b.WriteConfiguration(input.writefile, iteration);
    }
  
  output.close();

  std::cout << "b.pf = " << b.pf << std::endl;
  std::cout << "b.pressure = " << b.pressure << std::endl;
  std::cout << "b.collisionrate = " << b.collisionrate << std::endl;
    
  return 0;
}
