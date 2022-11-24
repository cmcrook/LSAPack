#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <chrono>

#include "read_input.hpp"
#include <string.h>

#ifdef WIN32
#define strcasecmp _stricmp
#endif 

void load_dist_file(const char* wconfigfile, std::vector<double>& sizes, std::vector<int>& count) {
	std::cout << "Reading void size distribution from file..." << std::endl;

	//Start reading input file
	std::ifstream file;
	file.open(wconfigfile, std::ifstream::in);
	if (!file.is_open()) {
		throw;
	}

	int line_num = 0;
	std::string tp;
	std::cout << "\t" << "Sizes" << "           " << "Frequency" << std::endl;
	while (std::getline(file, tp)) {
		//Pad tp with extra space
		tp += " ";

		size_t pos = 0;
		int term_count = 0;
		while ((pos = tp.find_first_of(" \n\0")) != std::string::npos) {
			std::string term = tp.substr(0, pos);
			if (!term.empty()) {
				double val = std::stod(term);
				switch (term_count) {
				case 0:
					sizes.push_back(val);
					break;
				case 1:
					count.push_back(val);
					break;
				default:
					std::cout << "Too many terms in distribution on line " << line_num << std::endl;
					throw;
				}

				term_count++;
			}
			tp.erase(0, pos + 1);
		}

		std::cout << "\t" << std::left << std::setw(15) << std::setprecision(7) << sizes[line_num] << " " << count[line_num] << std::endl;

		line_num++;
	}

	file.close();
}


// Source File for input
int read_input::read(int argc, char* argv[])
{
	int error = 0;
	if (argc != 2) {
		std::cout << "Error, no input script provided!" << std::endl;
		exit(-1);
	}
	else
	{
		std::ifstream infile;
		infile.open(argv[1]);
		if (!infile) {
			std::cout << "Can't open the input file " << argv[1] << std::endl;
			exit(-1);
		}
		else {
			std::cout << "Reading input file " << argv[1] << std::endl;
		}
		char buf[100], c;
		infile.get(buf, 100, '='); infile.get(c); infile >> dim;
		infile.get(buf, 100, '='); infile.get(c); infile >> seed;
		infile.get(buf, 100, '='); infile.get(c); infile >> eventspercycle;
		infile.get(buf, 100, '='); infile.get(c); infile >> initialpf;
		infile.get(buf, 100, '='); infile.get(c); infile >> maxpf;
		infile.get(buf, 100, '='); infile.get(c); infile >> temp;
		infile.get(buf, 100, '='); infile.get(c); infile >> growthrate;
		infile.get(buf, 100, '='); infile.get(c); infile >> maxpressure;
		infile.get(buf, 100, '='); infile.get(c); infile >> maxcollisionrate;
		infile.get(buf, 100, '='); infile.get(c); infile >> maxSizeChange;
		infile.get(buf, 100, '='); infile.get(c); infile >> citer;

		//Need to rewrite to read in list of sizes, fractions and masses from input file
		infile.get(buf, 100, '='); infile.get(c); infile >> hardwalls;
		infile.get(buf, 100, '='); infile.get(c);
		infile.width(NAME_LEN - 1); infile >> readfile;
		infile.get(buf, 100, '='); infile.get(c);
		infile.width(NAME_LEN - 1); infile >> distfile;
		infile.get(buf, 100, '='); infile.get(c);
		infile.width(NAME_LEN - 1); infile >> writefile;
		infile.get(buf, 100, '='); infile.get(c);
		infile.width(NAME_LEN - 1); infile >> datafile;

		if (infile.eof())
		{
			std::cout << "Error reading input file " << argv[1] << std::endl;
			exit(-1);
		}

		if (seed == 0) {
			seed = std::chrono::system_clock::now().time_since_epoch().count();
			std::cout << "No seed provided, randomly chose " << static_cast<int>(seed) << std::endl;
		}

		if (strcasecmp(distfile, "none") != 0) {
			load_dist_file(distfile, particle_sizes, particle_counts);

			N = 0;
			for (auto c : particle_counts) {
				N += c;
			}
		}
		else {
			std::cout << "No distribute file specified!" << std::endl;
			throw;
		}

		std::cout << "\tdim : " << dim << std::endl;
		std::cout << "\tseed : " << seed << std::endl;
		std::cout << "\teventspercycle : " << eventspercycle << std::endl;
		std::cout << "\tN : " << N << std::endl;
		std::cout << "\tinitialpf : " << initialpf << std::endl;
		std::cout << "\tmaxpf : " << maxpf << std::endl;
		std::cout << "\ttemp : " << temp << std::endl;
		std::cout << "\tgrowthrate : " << growthrate << std::endl;
		std::cout << "\tmaxpressure : " << maxpressure << std::endl;
		std::cout << "\tmaxcollisionrate : " << maxcollisionrate << std::endl;
		std::cout << "\tmaxSizeChange: " << maxSizeChange << std::endl;
		std::cout << "\tciter: " << citer << std::endl;
		std::cout << "\thardwallBC : " << hardwalls << std::endl;
		std::cout << "\treadfile : " << readfile << std::endl;
		std::cout << "\tdistfile : " << distfile << std::endl;
		std::cout << "\twritefile : " << writefile << std::endl;
		std::cout << "\tdatafile : " << datafile << std::endl;
	}

	return error;
}
