#pragma once
/*
   Packing of hard spheres via molecular dynamics
   Developed by Monica Skoge, 2006, Princeton University
   Contact: Aleksandar Donev (adonev@math.princeton.edu) with questions
   This code may be used, modified and distributed freely.
   Please cite:

   "Packing Hyperspheres in High-Dimensional Euclidean Spaces"
	M. Skoge, A. Donev, F. H. Stillinger and S. Torquato, 2006

   if you use these codes.
*/

//-----------------------------------------------------------------------------
// Box maker
//---------------------------------------------------------------------------

#include <vector>
#include <math.h>

#include "vec.hpp"
#include "grid_field.hpp"
#include "event.hpp"
#include "sphere.hpp"
#include "heap.hpp"


#define PI     3.141592653589793238462643
#define SIZE   1.0            // size of box
#define VOLUMESPHERE pow(PI,((double)(DIM))/2.)/exp(lgamma(1+((double)(DIM))/2.)) // volume prefactor for sphere
#define DBL_EPSILON  2.2204460492503131e-016 // smallest # such that 1.0+DBL_EPSILON!=1.0
#define M 1.0

//---------------------------------------------------------------------------
// Class neighbor
//---------------------------------------------------------------------------
class neighbor
{
public:
	int i;

	neighbor(int i_i);

public:
	virtual void Operation(int j, vector<DIM, int>& pboffset) = 0;
};


class Box {
public:
	// constructor and destructor
	Box(int N_i, double r_i, double growthrate_i, double maxpf_i,
		std::vector<double> bidispersityratio, std::vector<double> bidispersityfraction,
		std::vector<double> massratio, int hardwallBC);
	~Box();

	// Creating configurations
	int optimalngrids();
	void createSpheres(double temp);
	void createSphere(int Ncurrent, double radius, double growth_rate, double mass, int species);
	double velocity(double temp);
	void thermalize(double temp);
	void SetInitialEvents();
	void recreateSpheres(const char* filename, double temp);
	void readPositions(const char* filename);
	void assignCells();

	// Predicting next event
	Event findNextEvent(int i);
	void collisionChecker(Event c);
	Event findNextTransfer(int i);
	Event findNextCollision(int i);
	void forAllNeighbors(int, vector<DIM, int>, vector<DIM, int>, neighbor&);
	void predictCollision(int i, int j, vector<DIM, int> pboffset, double& ctime, int& cpartner, vector<DIM, int>& cpartnerpboffset);
	double calculateCollision(int i, int j, vector<DIM>  pboffset);
	double quadraticFormula(double a, double b, double c);

	// Processing an event
	void process(int n);
	void processEvent();
	void collision(Event e);
	void transfer(Event e);
	void updateCell(int i, vector<DIM, int>& celli);
	void synchronize(bool rescale);
	void changeNgrids(int newngrids);

	// Debugging
	void trackPositions();
	void outputEvents();
	void outputCells();
	void GetInfo();
	int checkSphereDiameters();

	// Statistics
	double Energy();
	double packingFraction();
	void printStatistics();
	void runTime();
	void writeConfiguration(const char* wconfigfile, int iteration);


	//variables

	const int N;                   // number of spheres

	int ngrids;                    // number of cells in one direction
	double maxpf;
	double growthrate;             // growth rate of the spheres
	double r;                      // radius, defined at gtime = 0
	double gtime;                  // this is global clock
	double rtime;                  // reset time, total time = rtime + gtime
	double collisionrate;          // average rate of collision between spheres
	double maxSizeChange;
	std::vector<double> particle_sizes;      // ratio of sphere radii
	std::vector<double> particle_fraction;   // fraction of smaller spheres
	std::vector<double> particle_mass;              // ratio of sphere masses
	std::vector<int> particle_count;					//count of each particle size
	int hardwallBC;                // =0 for periodic BC, =1 for hard wall

	// statistics
	double pressure;               // pressure
	double xmomentum;              // exchanged momentum
	double pf;                     // packing fraction
	double energy;                 // kinetic energy
	double energychange;
	int ncollisions;               // number of collisions
	int ntransfers;                // number of transfers
	int nchecks;                   // number of checks
	int neventstot;                // total number of events 
	int ncycles;                   // counts # cycles for output

	time_t start, error, end;      // run time of program

	// arrays
	Sphere* s;                      // array of spheres
	grid_field<DIM, int> cells; // array that keeps track of spheres in each cell
	int* binlist;                   // linked-list for cells array
	heap heap;                         // event heap
	vector<DIM>* x;                 // positions of spheres.used for graphics
};


//---------------------------------------------------------------------------
// Predicts collisions, inherits neighbor operation
//---------------------------------------------------------------------------
class Collision : public neighbor {
public:

	Box* b;
	double ctime;
	int cpartner;
	vector<DIM, int> cpartnerpboffset;

public:
	Collision(int i_i, Box* b);

	virtual void Operation(int j, vector<DIM, int>& pboffset);
};