//  Molecular dynamics simulation of hardspheres

#define strncasecmp _strnicmp
#define strcasecmp _stricmp

#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <time.h>
#include <string.h>
#include <iomanip>

#include "box.hpp"
#include "sphere.hpp"
#include "event.hpp"
#include "heap.hpp"
#include "read_input.hpp"


int main(int argc, char** argv)
{
	read_input input;
	int error = input.read(argc, argv);
	if (error) return error;

	double d, r;   // initial diameter and radius of spheres

	if (strcasecmp(input.readfile, "new") == 0)
		input.readfile[0] = 0;

	if (input.readfile[0]) // read in existing configuration
	{
		// read the header
		std::ifstream infile(input.readfile);
		if (!infile)
		{
			std::cout << "error, can't open " << input.readfile << std::endl;
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
			r = d / 2.;
			std::cout << "r = " << r << std::endl;
		}
	}
	else // create a new configuration
	{
		r = pow(input.initialpf * pow(SIZE, DIM) / (input.N * VOLUMESPHERE), 1.0 / ((double)(DIM)));
	}

	//Temporary radii distribution for testing
	std::vector<double> particle_masses(input.particle_sizes.size(), 1.0);

	Box b(input.N, r, input.growthrate, input.maxpf, input.particle_sizes, input.particle_fractions, particle_masses, input.hardwallBC);

	std::cout << "ngrids = " << b.ngrids << std::endl;
	std::cout << "DIM = " << DIM << std::endl;

	if (input.readfile[0])
	{
		std::cout << "Reading in positions of spheres" << std::endl;
		b.recreateSpheres(input.readfile, input.temp);
	}
	else
	{
		std::cout << "Creating new positions of spheres" << std::endl;
		b.createSpheres(input.temp);
	}

	std::ofstream output(input.datafile);
	output.precision(16);
	output.setf(std::ios::fixed, std::ios::floatfield);
	
	std::cout << std::setw(10) << "Running LSA..." << std::endl;
	std::cout << "\t" << std::setw(8) << "Iter" << " " << std::setw(15) << "Collision Rate" << " " << std::setw(10) << "Pf" << " " << std::setw(15) << "Press" << " " << std::setw(10) << "Max dR" << std::endl;

	int iteration = 0;
	b.writeConfiguration(input.writefile, iteration);
	while ((b.pf < input.maxpf) && (b.maxSizeChange > input.maxSizeChange))	{
		b.process(input.eventspercycle * input.N);
		output << b.pf << " " << b.pressure << " " << b.collisionrate << " " << b.neventstot << " " << std::endl;

		b.synchronize(true);
		std::cout << "\t" << std::left << std::setw(8) << iteration << " " << std::left << std::setw(15) << std::setprecision(7) << std::scientific << b.collisionrate << " " << std::left << std::setw(10) << std::fixed << b.pf << "  " << std::left << std::setw(15) << std::scientific << b.pressure << " " << std::left << std::setw(10) << b.maxSizeChange << std::endl;

		iteration++;
		b.writeConfiguration(input.writefile, iteration);
	}

	/*b.writeConfiguration(input.writefile, iteration);*/

	output.close();

	b.printStatistics();

	return 0;
}
