//  Molecular dynamics simulation of hardspheres
//  Modified by Cameron Crook

#ifdef WIN32
#define strcasecmp _stricmp
#endif 

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <fstream>
#include <math.h>
#include <time.h>
#include <string.h>

#include "box.hpp"
#include "sphere.hpp"
#include "event.hpp"
#include "heap.hpp"
#include "read_input.hpp"
#include "recorder.hpp"

//Compute RDF
void pack2rdf(Box& box, int num_voids, int num_bins, double cutoff) {

	//Determine image offsets
	std::vector<double> xOff, yOff, zOff;
	xOff.push_back(0.0);
	yOff.push_back(0.0);
	zOff.push_back(0.0);

	if (!box.hardwallBC) {
		xOff.push_back(-1.0);
		xOff.push_back(1.0);
		yOff.push_back(-1.0);
		yOff.push_back(1.0);
		zOff.push_back(-1.0);
		zOff.push_back(1.0);
	}

	double bin_width = cutoff / num_bins;
	std::vector<int> bins(num_bins, 0);

	//Currently does not include periodicity
	for (int i = 0; i < num_voids; i++) {
		for (int j = i + 1; j < num_voids; j++) {
			vector<DIM> sep;
			for(int k = 0; k < DIM; k++)
				sep[k] = box.s[i].x[k] - box.s[j].x[k];

			for (const auto& xx : xOff) {
				for (const auto& yy : yOff) {
					for (const auto& zz : zOff) {
						//FG::Vector3d image = sep + FG::Vector3d(aabb.getWidth(), aabb.getHeight(), aabb.getLength()) * FG::Vector3d(xx, yy, zz);
						vector<DIM> image;
						image[0] = SIZE * xx  + sep[0];
						image[1] = SIZE * yy + sep[1];
						image[2] = SIZE * zz + sep[2];

						double sqrDist = image[0] * image[0] + image[1] * image[1] + image[2] * image[2];
						if (sqrDist < cutoff * cutoff) {
							int bin_i = floor(std::sqrt(sqrDist) / bin_width);
							bins[bin_i]++;
						}
					}
				}
			}
		}
	}

	//Format to string
	double c = 2.0 * SIZE * SIZE * SIZE / (num_voids * num_voids * 4.0 * PI * bin_width);
	std::ostringstream ss;
	for (int i = 0; i < num_bins; i++) {
		double r = (i + 0.5) * bin_width;
		ss << c * bins[i] / (r * r) << std::endl;
	}

	//Write rdf to file
	std::ofstream file;
	file.open("PACK_RDF");
	file << ss.str();
	file.close();
}

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

	Box b(input.N, r, input.growthrate, input.maxpf, input.particle_sizes, input.particle_fractions, particle_masses, input.hardwallBC, input.seed);

	std::cout << "ngrids = " << b.ngrids << std::endl;
	std::cout << "DIM = " << DIM << std::endl;

	if (input.readfile[0])
	{
		std::cout << "Reading in sphere positions from file" << std::endl;
		b.recreateSpheres(input.readfile, input.temp);
	}
	else
	{
		std::cout << "Creating new sphere  positions" << std::endl;
		b.createSpheres(input.temp);
	}

	std::ofstream output(input.datafile);
	output.precision(16);
	output.setf(std::ios::fixed, std::ios::floatfield);

	std::cout << std::setw(10) << "Running LSA..." << std::endl;
	std::stringstream stats;
	stats << "\t" << std::left << std::setw(5) << "Iter" 
		  << " " << std::left << std::setw(10) << "Pf" 
		  << " " << std::left << std::setw(15) << "Collision Rate" 
		  << " " << std::left << std::setw(15) << "Press" 
		  << " " << std::left << std::setw(15) << "Max dR" 
		  << " " << std::left << std::setw(7) << "Ngrids"
		  << " " << std::left << std::setw(10) << "Sim Time" << std::endl;
	output << stats.str();
	std::cout << stats.str();

	int iteration = 0;
	//b.writeLAMMPSDump("pack.dump", iteration);
	b.writeConfiguration(input.writefile, iteration);
	while (b.pf < input.maxpf)	{
		iteration++; //0 is initial configuration
		
		{
			PROFILE_SCOPE("Process");
			b.process(input.eventspercycle * input.N);
		}

		{
			PROFILE_SCOPE("Synchronize");
			b.synchronize(true);
		}

		{
			PROFILE_SCOPE("Write configuration");
			stats.str("");
			stats << "\t" << std::left << std::setw(5) << iteration 
				<< " " << std::left << std::setw(10) << std::fixed << b.pf 
				<< " " << std::left << std::setw(15) << std::setprecision(7) << std::scientific << b.collisionrate 
				<< " " << std::left << std::setw(15) << std::scientific << b.pressure 
				<< " " << std::left << std::setw(15) << b.maxSizeChange 
				<< " " << std::left << std::setw(7) << b.ngrids
				<< " " << std::left << std::setw(10) << b.rtime << std::endl;
			output << stats.str();
			std::cout << stats.str();

			b.writeConfiguration(input.writefile, iteration);
		}

		if (b.maxSizeChange < input.maxSizeChange)
			b.citer++;
		else
			b.citer = 0;

		if (b.citer >= input.citer)
			break;
	}

	//Output is now a custom lammps dump file, rdf can be calculated in ovito
	//pack2rdf(b, input.N, 100, SIZE);

	b.writeLAMMPSDump("pack.dump", iteration);

	output.close();

	b.printStatistics();

	Recorder::printStats();

	return 0;
}
