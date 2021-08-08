/*
   Updated July 24, 2009 to include hardwall boundary condition option, as
   well as polydisperse spheres.

   Packing of hard spheres via molecular dynamics
   Developed by Monica Skoge, 2006, Princeton University
   Contact: Monica Skoge (mskoge.princeton.edu) with questions
   This code may be used, modified and distributed freely.
   Please cite:

   "Packing Hyperspheres in High-Dimensional Euclidean Spaces"
	M. Skoge, A. Donev, F. H. Stillinger and S. Torquato, 2006

   if you use these codes.
*/

#include "box.hpp"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iomanip>
#include <vector>
#include <algorithm>

//==============================================================
//==============================================================
//  Class Box: Fills box with hardspheres to given packing fraction
//  and evolves spheres using molecular dynamics!
//==============================================================
//==============================================================

//==============================================================
// Constructor
//==============================================================
Box::Box(int N_i, double r_i, double growthrate_i, double maxpf_i,
	 std::vector<double> bidispersityratio_i, std::vector<double> bidispersityfraction_i,
	std::vector<double> massratio_i, int hardwallBC_i, int seed) :
	r(r_i),
	N(N_i),
	growthrate(growthrate_i),
	heap(N_i + 1),
	maxpf(maxpf_i),
	particle_sizes(bidispersityratio_i),
	particle_fraction(bidispersityfraction_i),
	particle_mass(massratio_i),
	hardwallBC(hardwallBC_i),
	seed(seed)
{
	gtime = 0.0;
	rtime = 0.0;
	ncollisions = 0;
	ntransfers = 0;
	nchecks = 0;
	neventstot = 0;
	ncycles = 0;
	xmomentum = 0.0;
	pressure = 0.0;
	collisionrate = 0.0;
	maxSizeChange = DBL_MAX;

	ngrids = optimalngrids();
	cells.set_size(ngrids);
	cells.initialize(-1);      // initialize cells to -1

	s = new Sphere[N];
	binlist = new int[N];
	x = new vector<DIM>[N];
	heap.s = s;

	srand(seed);        // initialize the random number generator
	for (int i = 0; i < N; i++)    // initialize binlist to -1
		binlist[i] = -1;

	time(&start);
}


//==============================================================
// Destructor
//==============================================================
Box::~Box() {
	delete[] s;
	delete[] binlist;
	delete[] x;
}

// ReadFile
void Box::readPositions(const char* filename) {
	// open file to read in arrays
	std::ifstream infile(filename);

	infile.ignore(256, '\n');  // ignore the dim line
	infile.ignore(256, '\n');  // ignore the #sphere 1 line
	infile.ignore(256, '\n');  // ignore the #sphere line
	infile.ignore(256, '\n');  // ignore the diameter line
	infile.ignore(1000, '\n'); // ignore the 100 010 001 line
	infile.ignore(256, '\n');  // ignore the T T T line

	for (int i = 0; i < N; i++)
	{
		infile >> s[i].r;      // read in radius    
		infile >> s[i].gr;     // read in growth rate
		infile >> s[i].m;      // read in mass
		for (int k = 0; k < DIM; k++)
			infile >> s[i].x[k]; // read in position 
	}
	infile.close();
}

// Recreates all N spheres at random positions
void Box::recreateSpheres(const char* filename, double temp) {
	readPositions(filename);  // reads in positions of spheres
	thermalize(temp);      // gives spheres initial velocities
	assignCells();            // assigns spheres to cells
	setInitialEvents();
}

// Creates all N spheres at random positions
void Box::createSpheres(double temp) {
	particle_count.resize(particle_fraction.size(), 0);
	int species = 0;
	double frac = 0.0;
	frac += particle_fraction[species];
	double radius = particle_sizes[species];
	double growth_rate = growthrate * particle_sizes[species]; //Shouldn't the growthrate be the same for all particles
	double mass = particle_mass[species];

	int Ncurrent = 0;
	for (int i = 0; i < N; i++)
	{
		createSphere(Ncurrent, r * radius, growth_rate, mass, species);
		Ncurrent++;
		particle_count[species]++;

		if (Ncurrent >= N * frac && species < particle_sizes.size() - 1) {
			species++;
			frac += particle_fraction[species];
			radius = particle_sizes[species];
			growth_rate = growthrate * particle_sizes[species];
			mass = particle_mass[species];
		}
	}
	if (Ncurrent != N)
		std::cout << "problem! only made " << Ncurrent << " out of " << N << " desired spheres" << std::endl;

	thermalize(temp);
	setInitialEvents();
}


//==============================================================
// Creates a sphere of random radius r at a random unoccupied position 
//==============================================================
void Box::createSphere(int Ncurrent, double radius, double growth_rate, double mass, int species)
{
	int keeper;    // boolean variable: 1 means ok, 0 means sphere already there
	int counter = 0;   // counts how many times sphere already exists
	vector<DIM> xrand;  // random new position vector
	double d = 0.0;

	while (counter < 1000)
	{
		keeper = 1;

		for (int k = 0; k < DIM; k++)
			xrand[k] = ((double)rand() / (double)RAND_MAX)* SIZE;

		for (int i = 0; i < Ncurrent; i++) {
			d = 0.0;
			if (hardwallBC) {
				for (int k = 0; k < DIM; k++) {
					d += (xrand[k] - s[i].x[k]) * (xrand[k] - s[i].x[k]);
				}
			} else {
				for (int k = 0; k < DIM; k++) {
					if ((xrand[k] - s[i].x[k]) * (xrand[k] - s[i].x[k]) > SIZE * SIZE / 4.) {
						if (xrand[k] > SIZE / 2.)  // look at right image
							d += (xrand[k] - (s[i].x[k] + SIZE)) * (xrand[k] - (s[i].x[k] + SIZE));
						else                     // look at left image
							d += (xrand[k] - (s[i].x[k] - SIZE)) * (xrand[k] - (s[i].x[k] - SIZE));
					} else
						d += (xrand[k] - s[i].x[k]) * (xrand[k] - s[i].x[k]);
				}
			}

			// overlapping!
			if (d <= (radius + s[i].r) * (radius + s[i].r)) {
				keeper = 0;
				counter++;
				break;
			}
		}

		if (hardwallBC)
		{
			for (int k = 0; k < DIM; k++)    // check if overlapping wall
			{
				if ((xrand[k] <= radius) || (SIZE - xrand[k] <= radius)) // touching wall
				{
					keeper = 0;
					counter++;
					break;
				}
			}
		}

		if (keeper == 1)
			break;

		if (counter >= 1000) {
			std::cout << "counter >= 1000" << std::endl;
			exit(-1);
		}
	}

	//Edge case coordinate lies on right edge of box
	for (int k = 0; k < DIM; k++) {
		if (xrand.x[k] == SIZE) {
			//std::cout << xrand.x[k] << std::endl;
			xrand.x[k] -= SIZE;
		}
	}
	// now convert xrand into index vector for cells
	vector<DIM, int> cell;
	cell = vector<DIM>::integer(xrand * ((double)(ngrids)) / SIZE);

	s[Ncurrent] = Sphere(Ncurrent, xrand, cell, gtime, radius, growth_rate, mass, species);

		//first check to see if entry at cell
		/*std::cout << xrand[0] << " " << xrand[1] << " " << xrand[2] << std::endl;
		std::cout << cell[0] << " " << cell[1] << " " << cell[2] << std::endl << std::endl;*/
		
		if (cells.get(cell) == -1) //if yes, add Ncurrent to cells gridfield
			cells.get(cell) = Ncurrent;

		else  // if no, add i to right place in binlist
		{
			int iterater = cells.get(cell); // now iterate through to end and add Ncurrent
			int pointer = iterater;
			while (iterater != -1)
			{
				pointer = iterater;
				iterater = binlist[iterater];
			}
			binlist[pointer] = Ncurrent;
		}
}


//==============================================================
// Assign cells to spheres read in from existing configuration
//==============================================================
void Box::assignCells() {
	for (int i = 0; i < N; i++)
	{
		// now convert x into index vector for cells
		vector<DIM, int> cell;
		cell = vector<DIM>::integer(s[i].x * ((double)(ngrids)) / SIZE);
		s[i].cell = cell;

		//first check to see if entry at cell
		//std::cout << s[i].x[0] << " " << s[i].x[1] << " " << s[i].x[2] << std::endl;
		//std::cout << cell[0] << " " << cell[1] << " " << cell[2] << std::endl << std::endl;
		if (cells.get(cell) == -1) //if yes, add Ncurrent to cells gridfield
			cells.get(cell) = i;

		else  // if no, add i to right place in binlist
		{
			int iterater = cells.get(cell); // now iterate through to end and add Ncurrent
			int pointer = iterater;
			while (iterater != -1)
			{
				pointer = iterater;
				iterater = binlist[iterater];
			}
			binlist[pointer] = i;
		}
	}
}

// Assigns initial velocities from Max/Boltz dist.
void Box::thermalize(double temp) {
	for (int i = 0; i < N; i++)
	{
		for (int k = 0; k < DIM; k++)
		{
			if (temp == 0.0)
				s[i].v[k] = 0.0;
			else
				s[i].v[k] = velocity(temp);
		}
	}
}

// Velocity, gives a single velocity from Max/Boltz dist.
double Box::velocity(double temp) {
	double random;                       // random number between -0.5 and 0.5
	double sigmasquared = temp;    // Assumes M = mass of sphere = 1
	double sigma = sqrt(sigmasquared); // variance of Gaussian
	double stepsize = 1000.;           // stepsize for discretization of integral
	double vel = 0.0;                  // velocity
	double dv = sigma / stepsize;
	double p = 0.0;

	random = (double)rand() / (double)RAND_MAX - 0.5;
	if (rand < 0)
	{
		random = -random;
		dv = -dv;
	}

	while (fabs(p) < random) // integrate until the integral equals rand
	{
		p += dv * 0.39894228 * exp(-vel * vel / (2. * sigmasquared)) / sigma;
		vel += dv;
	}
	return vel;
}


//==============================================================
// Finds next events for all spheres..do this once at beginning
//==============================================================
void Box::setInitialEvents()
{
	for (int i = 0; i < N; i++)  // set all events to checks
	{
		Event e(gtime, i, INF);
		s[i].nextevent = e;
		heap.insert(i);
	}
}

// Finds next event for sphere i 
Event Box::findNextEvent(int i)
{
	double outgrowtime;
	outgrowtime = (0.5 / ngrids - (s[i].r + s[i].gr * gtime)) / s[i].gr + gtime;

	Event t = findNextTransfer(i);
	Event c = findNextCollision(i);

	//std::cout << outgrowtime << " " << c.time << " " << t.time << std::endl;
	//std::cout << (outgrowtime < c.time) << " " << (outgrowtime < t.time) << " " << (ngrids > 1) << std::endl;

	//Temporary fix attempt to prevent outgrowtime from being negative
	outgrowtime = std::max(outgrowtime, 0.0);

	if ((outgrowtime < c.time) && (outgrowtime < t.time) && (ngrids > 1))
	{
		Event o = Event(outgrowtime, i, INF - 1);
		return o;
	}

	if ((c.time < t.time) && (c.j == INF)) // next event is check at DBL infinity
		return c;
	else if (c.time < t.time) // next event is collision!
	{
		collisionChecker(c);
		return c;
	}
	else // next event is transfer!
		return t;
}


//==============================================================
// Checks events of predicted collision partner to keep collisions
// symmetric
//==============================================================
void Box::collisionChecker(Event c)
{
	int i = c.i;
	int j = c.j;
	Event cj(c.time, j, i, c.v * (-1));

	// j should have NO event before collision with i!
	if (!(c.time < s[j].nextevent.time))
		std::cout << i << " " << j << " error collchecker, s[j].nextevent.time= " << s[j].nextevent.time << " " << s[j].nextevent.j << ", c.time= " << c.time << std::endl;

	int k = s[j].nextevent.j;
	if ((k < N) && (k != i)) // j's next event was collision so give k a check
		s[k].nextevent.j = INF;

	// give collision cj to j
	s[j].nextevent = cj;
	heap.upheap(heap.index[j]);
}


//==============================================================
// Find next transfer for sphere i 
//==============================================================
Event Box::findNextTransfer(int i)
{
	double ttime = dblINF;
	int wallindex = INF;   // -(k+1) left wall, (k+1) right wall

	vector<DIM> xi = s[i].x + s[i].v * (gtime - s[i].lutime);
	vector<DIM> vi = s[i].v;

	for (int k = 0; k < DIM; k++)
	{
		double newtime;
		if (vi[k] == 0.)
			newtime = dblINF;
		else if (vi[k] > 0)  // will hit right wall, need to check if last wall
		{
			if ((hardwallBC) && (s[i].cell[k] == ngrids - 1))
				newtime = ((double)(s[i].cell[k] + 1) * SIZE / ((double)(ngrids))
					   - (xi[k] + s[i].r + s[i].gr * gtime)) / (vi[k] + s[i].gr);
			else
				newtime = ((double)(s[i].cell[k] + 1) * SIZE / ((double)(ngrids))
					   - xi[k]) / (vi[k]);

			if (newtime < ttime)
			{
				wallindex = k + 1;
				ttime = newtime;
			}
		}
		else if (vi[k] < 0)  // will hit left wall
		{
			if ((hardwallBC) && (s[i].cell[k] == 0))
				newtime = ((double)(s[i].cell[k]) * SIZE / ((double)(ngrids))
					   - (xi[k] - (s[i].r + s[i].gr * gtime))) / (vi[k] - s[i].gr);
			else
				newtime = ((double)(s[i].cell[k]) * SIZE / ((double)(ngrids))
					   - xi[k]) / (vi[k]);

			if (newtime < ttime)
			{
				wallindex = -(k + 1);
				ttime = newtime;
			}
		}
	}


	if (ttime < 0)
		ttime = 0;
	// make the event and return it
	Event e = Event(ttime + gtime, i, wallindex + DIM + N + 1);
	return e;
}


//==============================================================
// Check all nearest neighbor cells for collision partners
//==============================================================
void Box::forAllNeighbors(int i, vector<DIM, int> vl, vector<DIM, int> vr,
			  Neighbor& operation)
{
	vector<DIM, int> cell = s[i].cell;

	// now iterate through nearest neighbors
	vector<DIM, int> offset;          // nonnegative neighbor offset
	vector<DIM, int> pboffset;        // nearest image offset

	vector<DIM, int> grid;

	int ii;

	grid = vl;
	while (1)
	{
		//if (vr[0] > 1)
		//std::cout << grid << "..." << cell+grid << "\n";
		for (int k = 0; k < DIM; k++)
		{
			offset[k] = grid[k] + ngrids;  // do this so no negatives 	
			if (cell[k] + grid[k] < 0) //out of bounds to left
				pboffset[k] = -1;
			else if (cell[k] + grid[k] >= ngrids) // out of bounds to right
				pboffset[k] = 1;
			else
				pboffset[k] = 0;
		}
		int j = cells.get((cell + offset) % ngrids);
		while (j != -1)
		{
			operation.Operation(j, pboffset);
			j = binlist[j];
		}

		// A. Donev:     
		// This code makes this loop dimension-independent
		// It is basically a flattened-out loop nest of depth DIM
		for (ii = 0; ii < DIM; ii++)
		{
			grid[ii] += 1;
			if (grid[ii] <= vr[ii]) break;
			grid[ii] = vl[ii];
		}
		if (ii >= DIM) break;
	}
}


//==============================================================
// PredictCollision
//==============================================================
void Box::predictCollision(int i, int j, vector<DIM, int> pboffset,
				 double& ctime, int& cpartner,
				 vector<DIM, int>& cpartnerpboffset)
{
	double ctimej;

	if (i != j)
	{
		ctimej = calculateCollision(i, j, pboffset.Double()) + gtime;

		if (ctimej < gtime)
			std::cout << "error in find collision ctimej < 0" << std::endl;

		if ((ctimej < ctime) && (ctimej < s[j].nextevent.time))
		{
			ctime = ctimej;
			cpartner = j;
			cpartnerpboffset = pboffset;
		}
	}
}


//==============================================================
// Find next collision
//==============================================================
Event Box::findNextCollision(int i) {
	Collision cc(i, this);

	vector<DIM, int> vl, vr;

	for (int k = 0; k < DIM; k++)  // check all nearest neighbors
	{
		vl[k] = -1;
		vr[k] = 1;
	}

	forAllNeighbors(i, vl, vr, cc);

	Event e;
	if (cc.cpartner == i)  // found no collisions in neighboring cells
	{
		if (cc.ctime != dblINF)
			std::cout << "ctime != dblINF" << std::endl;
		e = Event(dblINF, i, INF);  // give check at double INF
	}
	else
		e = Event(cc.ctime, i, cc.cpartner, cc.cpartnerpboffset);

	return e;
}

// Calculates collision time between i and image of j using quadratic formula
double Box::calculateCollision(int i, int j, vector<DIM> pboffset) {
	if ((hardwallBC) && (pboffset.norm_squared() > 1E-12))
	{
		return dblINF;
	}
	else
	{
		// calculate updated position and velocity of i and j
		vector<DIM> xi = s[i].x + s[i].v * (gtime - s[i].lutime);
		vector<DIM> vi = s[i].v;
		vector<DIM> xj = s[j].x + pboffset * SIZE + s[j].v * (gtime - s[j].lutime);
		vector<DIM> vj = s[j].v;

		double r_now_i = s[i].r + gtime * s[i].gr;
		double r_now_j = s[j].r + gtime * s[j].gr;

		double A, B, C;
		A = vector<DIM>::norm_squared(vi - vj) - (s[i].gr + s[j].gr) * (s[i].gr + s[j].gr);
		B = vector<DIM>::dot(xi - xj, vi - vj) - (r_now_i + r_now_j) * (s[i].gr + s[j].gr);
		C = vector<DIM>::norm_squared(xi - xj) - (r_now_i + r_now_j) * (r_now_i + r_now_j);

		if (C < -1E-12 * (r_now_i + r_now_j))
		{
			std::cout << "error, " << i << " and " << j <<
				" are overlapping at time " << gtime << std::endl;
			std::cout << "A, B, C = " << A << " " << " " << B <<
				" " << " " << C << std::endl;
			if (checkSphereDiameters() > 0)
				std::cout << "a sphere has grown greater than unit cell" <<
				std::endl;
			else
				std::cout << "unknown error" << std::endl;
			exit(-1);
		}

		return quadraticFormula(A, B, C);
	}
}

// Quadratic Formula ax^2 + bx + c = 0
double Box::quadraticFormula(double a, double b, double c)
{
	double x = dblINF;
	double xpos;
	double xneg;
	double det = b * b - a * c;

	if (c <= 0.)
	{
		if (b < 0.) // spheres already overlapping and approaching
		{
			//std::cout << "spheres overlapping and approaching" << std::endl;
			//std::cout << "# events= " << neventstot << std::endl;
			x = 0.;
		}
	}
	else if (det > -10. * DBL_EPSILON)
	{
		if (det < 0.)  // determinant can be very small for double roots
			det = 0.;
		if (b < 0.)
			x = c / (-b + sqrt(det));
		else if ((a < 0.) && (b > 0.))
			x = -(b + sqrt(det)) / a;
		else
			x = dblINF;
	}
	return x;
}

// Returns first event
void Box::processEvent() {
	neventstot++;
	// Extract first event from heap
	int i = heap.extractmax();
	Event e = s[i].nextevent; // current event
	Event f;                  // replacement event

	if ((e.j >= 0) && (e.j < N))  // collision!
	{
		ncollisions++;
		//std::cout << "collision between " << e.i << " and " << e.j << " at time " << e.time << std::endl;
		collision(e);
		f = findNextEvent(i);
		s[i].nextevent = f;
		heap.downheap(1);
		if (f.time < e.time)
		{
			std::cout << "error, replacing event with < time" << std::endl;
			exit(-1);
		}

		// make sure collision was symmetric and give j a check
		if ((s[e.j].nextevent.j != i) || (s[e.j].nextevent.time != gtime))
		{
			std::cout << "error collisions not symmetric" << std::endl;
			std::cout << "collision between " << e.i << " and " << e.j << " at time " << e.time << std::endl;
			std::cout << "but " << e.j << " thinks it has " << s[e.j].nextevent.j << " " << s[e.j].nextevent.time << std::endl;
			exit(-1);
		}
		else  // give j a check
			s[e.j].nextevent.j = INF;
	}
	else if (e.j == INF)      // check!  
	{
		nchecks++;
		//std::cout << "check for " << e.i << " at time " << e.time << std::endl;
		f = findNextEvent(i);
		s[i].nextevent = f;
		heap.downheap(1);
	}
	else if (e.j == INF - 1)      // sphere outgrowing unit cell, decrease ngrids!
	{
		gtime = e.time;
		synchronize(false);
		ngrids -= 1;
		//ngrids = Optimalngrids();
		//std::cout << "need to reduce ngrids to " << ngrids << std::endl;
		changeNgrids(ngrids);
		heap.downheap(1);
	}
	else                    // transfer!
	{
		ntransfers++;
		//std::cout << "transfer for " << e.i << " at time " << e.time << std::endl;
		transfer(e);
		f = findNextEvent(i);
		s[i].nextevent = f;
		heap.downheap(1);
		//r = FindNextEvent(i, e.j-N-DIM-1);
		//if (f.time <= e.time)
		if (f.time < e.time)
		{
			std::cout << "error after transfer, replacing new event with < time" << " " << std::endl;
			std::cout << std::setprecision(16) << "e.time= " << e.time << ", f.time= " << f.time << ", f.i= " << f.i << ", f.j= " << f.j << "e.i= " << e.i << ", e.j= " << e.j << std::endl;
			std::cout << std::setprecision(16) << "difference= " << e.time - f.time << std::endl;
			exit(-1);
		}
	}
}


// Processes a collision
void Box::collision(Event e)
{
	double ctime = e.time;
	int i = e.i;
	int j = e.j;
	vector<DIM, int> v = e.v;  // virtual image
	gtime = ctime;

	// Update positions and cells of i and j to ctime
	s[i].x += s[i].v * (gtime - s[i].lutime);
	s[j].x += s[j].v * (gtime - s[j].lutime);

	// Check to see if a diameter apart
	double r_sum = s[i].r + s[j].r + (s[i].gr + s[j].gr) * gtime;
	double distance = vector<DIM>::norm_squared(s[i].x - s[j].x - v.Double() * SIZE) - r_sum * r_sum;
	if (distance * distance > 10. * DBL_EPSILON)
		std::cout << "overlap " << distance << std::endl;

	s[i].lutime = gtime;
	s[j].lutime = gtime;

	vector<DIM, double> vipar;          // parallel comp. vi
	vector<DIM, double> vjpar;          // parallel comp. vj
	vector<DIM, double> viperp;         // perpendicular comp. vi
	vector<DIM, double> vjperp;         // perpendicular comp. vj

	// make unit vector out of displacement vector
	vector<DIM, double> dhat;
	dhat = s[i].x - s[j].x - v.Double() * SIZE;  // using image of j!!
	double dhatmagnitude = sqrt(dhat.norm_squared());
	dhat /= dhatmagnitude;

	vipar = dhat * vector<DIM>::dot(s[i].v, dhat);
	vjpar = dhat * vector<DIM>::dot(s[j].v, dhat);
	viperp = s[i].v - vipar;
	vjperp = s[j].v - vjpar;

	s[i].v = vjpar + dhat * (s[i].gr + s[j].gr) * 2 + viperp;
	s[j].v = vipar - dhat * (s[i].gr + s[j].gr) * 2 + vjperp;

	// momentum exchange
	double xvelocity;   // exchanged velocity
	xvelocity = vector<DIM>::dot(s[i].v - s[j].v, dhat) - (s[i].gr + s[j].gr);
	xmomentum += xvelocity * dhatmagnitude * s[i].m * s[j].m * 2 / (s[i].m + s[j].m);
}

// Transfer, takes care of boundary events too
void Box::transfer(Event e)
{
	gtime = e.time;
	int i = e.i;
	int j = e.j;
	int k = 0;           // dimension perpendicular to wall it crosses

	// update position and lutime (velocity doesn't change)
	s[i].x += s[i].v * (gtime - s[i].lutime);
	s[i].lutime = gtime;

	vector<DIM, int> celli;  // new cell for i
	celli = s[i].cell;  // this is not redundant

	// update cell
	if (j > N + DIM + 1)  // right wall
	{
		k = j - N - DIM - 2;
		celli[k] = s[i].cell[k] + 1;

		if (hardwallBC)
		{
			// if in right-most cell, reflect back
			if (s[i].cell[k] == ngrids - 1)
				s[i].v[k] = -s[i].v[k] - 4. * s[i].gr;
			else
			{
				if (ngrids > 1)
					updateCell(i, celli);
			}
		}
		else
		{
			// if in right-most cell, translate x and cell
			if (s[i].cell[k] == ngrids - 1)
			{
				s[i].x[k] -= SIZE;
				celli[k] -= ngrids;
			}
			if (ngrids > 1)
				updateCell(i, celli);
		}
	}
	else if (j < N + DIM + 1)  // left wall
	{
		k = -j + N + DIM;
		celli[k] = s[i].cell[k] - 1;

		if (hardwallBC)
		{
			// if in left-most cell, reflect back
			if (s[i].cell[k] == 0)
				s[i].v[k] = -s[i].v[k] + 4. * s[i].gr;
			else
			{
				if (ngrids > 1)
					updateCell(i, celli);
			}
		}
		else
		{
			// if in left-most cell, translate x and cell
			if (s[i].cell[k] == 0)
			{
				s[i].x[k] += SIZE;
				celli[k] += ngrids;
			}
			if (ngrids > 1)
				updateCell(i, celli);
		}
	}
	else
		std::cout << "error in Transfer" << std::endl;

}

// Updates cell of a sphere to time
void Box::updateCell(int i, vector<DIM, int>& celli) {
	if (celli == s[i].cell)
		std::cout << "error in update cell..shouldn't be the same" << std::endl;

	// delete i from cell array at cell

	if (cells.get(s[i].cell) == i)
	{
		if (binlist[i] == -1)
			cells.get(s[i].cell) = -1;
		else
		{
			cells.get(s[i].cell) = binlist[i];
			binlist[i] = -1;
		}
	}

	else if (cells.get(s[i].cell) == -1)
	{
		std::cout << "error " << i << " not in claimed cell UpdateCell" << std::endl;
		outputCells();
	}

	else  // if no, find i in binlist
	{
		int iterater = cells.get(s[i].cell);
		int pointer = iterater;
		while ((iterater != i) && (iterater != -1))
		{
			pointer = iterater;
			iterater = binlist[iterater];
		}
		if (iterater == -1)  // got to end of list without finding i
		{
			std::cout << "problem " << i << " wasn't in claimed, cell iterater = -1" << std::endl;
			outputCells();
		}
		else  // we found i!
		{
			binlist[pointer] = binlist[i];
			binlist[i] = -1;
		}
	}

	// now add i to cell array at celli
	s[i].cell = celli;

	//first check to see if entry at celli
	if (cells.get(celli) == -1) //if yes, add i to cells gridfield
		cells.get(celli) = i;
	else  // if no, add i to right place in binlist
	{
		int iterater = cells.get(celli);  // now iterate through to end and add i
		int pointer = iterater;
		while (iterater != -1)  // find the end of the list
		{
			pointer = iterater;
			iterater = binlist[iterater];
		}
		binlist[pointer] = i;
		binlist[i] = -1; // redundant
	}
}


//==============================================================
// Output event heap...purely used for debugging
//==============================================================
void Box::outputEvents()
{
	heap.print();
}


//==============================================================
// Output positions of spheres and their cells...purely used for debugging
//==============================================================
void Box::outputCells()
{
	for (int i = 0; i < N; i++)
		std::cout << i << " " << s[i].x << " " << s[i].v << " " << s[i].cell << std::endl;
}


//==============================================================
// Update positions...purely for graphical display
//==============================================================
void Box::trackPositions()
{
	for (int i = 0; i < N; i++)
		x[i] = s[i].x + s[i].v * (gtime - s[i].lutime);
}


//==============================================================
// Computes the total energy
//==============================================================
double Box::Energy() {
	double E = 0;
	for (int i = 0; i < N; i++)
		E += 0.5 * s[i].m * s[i].v.norm_squared();

	return E / N;
}


//==============================================================
// Calculates the packing fraction
//==============================================================
double Box::packingFraction() {
	double rfactor = 0.;
	for (int i = 0; i < N; i++)
	{
		rfactor += pow(s[i].r + gtime * s[i].gr, DIM);
	}
	double v = (rfactor * pow(sqrt(PI), DIM)) / (exp(lgamma(1. + ((double)(DIM)) / 2.)));
	return v / (pow(SIZE, DIM));
}


//==============================================================
// Checks to make sure all sphere diameters are less than dimension
// of unit cell
//==============================================================
int Box::checkSphereDiameters() {
	int offender = 0;
	for (int i = 0; i < N; i++)
	{
		if (s[i].r * 2 > 1. / ngrids) {
			offender = i;
			break;
		}
	}
	return offender;
}


//==============================================================
// Change ngrids
//==============================================================
void Box::changeNgrids(int newngrids) {
	cells.set_size(newngrids);
	cells.initialize(-1);      // initialize cells to -1
	for (int i = 0; i < N; i++)    // initialize binlist to -1
		binlist[i] = -1;
	assignCells();
	for (int i = 0; i < N; i++)
		s[i].nextevent = Event(0., i, INF);
	process(N);
}


//==============================================================
// Calculates the optimal ngrids for the initial configuration
// and assumes that ngrids gets updated (reduced) as the packing
// proceeds
//==============================================================
int Box::optimalngrids()
{
	// Calculates the optimal ngrids, assuming ngrids is not updated
	// automatically and is very conservative

	double largest_particle = 0.0;
	for (const auto& s : particle_sizes) {
		if (s > largest_particle) {
			largest_particle = s;
		}
	}

	double maxr = r * largest_particle + growthrate * gtime;

	//maxr = pow(exp(lgamma(1.+((double)(DIM))/2.))*maxpf/
	   //    (N*(bidispersityfraction + (1.-bidispersityfraction)*
		  // pow(bidispersityratio, DIM))), 
	   //    1./DIM)/sqrt(PI);

	int ngrid = (int)(1. / (2. * maxr));

	std::cout << "Max R: " << maxr << std::endl;
	std::cout << "Ngrids: " << ngrid << std::endl;
	return ngrid;
}


//==============================================================
// Processes n events
//==============================================================
void Box::process(int n)
{
	double deltat = gtime;
	for (int i = 0; i < n; i++)
	{
		processEvent();
		// std::cout << "Event: " << i << std::endl;
	}
	pf = packingFraction();   // packing fraction
	deltat = gtime - deltat;
	double oldenergy = energy;
	energy = Energy();        // kinetic energy

	energychange = ((oldenergy - energy) / oldenergy) * 100; // percent change in energy

	if (deltat != 0.)
	{
		pressure = 1 + xmomentum / (2. * energy * N * deltat);
		collisionrate = ((double)(ncollisions)) / deltat;
	}

	// reset to 0
	ncollisions = 0;
	ntransfers = 0;
	nchecks = 0;
	xmomentum = 0.;
	ncycles++;
}


//==============================================================
// Prints statistics for n events
//==============================================================
void Box::printStatistics() {
	std::cout << "seed = " << seed << std::endl;
	std::cout << "packing fraction = " << pf << std::endl;
	std::cout << "gtime = " << gtime << std::endl;
	std::cout << "total time = " << rtime + gtime << std::endl;
	std::cout << "kinetic energy = " << energy << std::endl;
	std::cout << "total # events = " << neventstot << std::endl;
	std::cout << "# events = " << ncollisions + ntransfers + nchecks << ", # collisions = " << ncollisions << ", # transfers = " << ntransfers << ", # checks =" << nchecks << std::endl;
	std::cout << "growthrate = " << growthrate << std::endl;
	std::cout << "collisionrate = " << collisionrate << std::endl;
	std::cout << "reduced pressure = " << pressure << std::endl;
	std::cout << "-----------------" << std::endl;
}


//==============================================================
// Updates spheres to gtime, synchronizes, and can change growth rate
//==============================================================
void Box::synchronize(bool rescale) {
	double vavg = sqrt(2. * M * energy);

	maxSizeChange = 0.0;
	double min_update_time = 0.0;
	for (int i = 0; i < N; i++)
	{
		if (min_update_time < (gtime - s[i].lutime)) {
			maxSizeChange = gtime - s[i].lutime;
		}

		s[i].x = s[i].x + s[i].v * (gtime - s[i].lutime);
		s[i].nextevent.time -= gtime;

		if (s[i].nextevent.time < 0.)
			std::cout << "error, event times negative after synchronization" << std::endl;
		if (rescale == true)   // give everyone checks
		{
			s[i].nextevent = Event(0., i, INF);
			s[i].v /= vavg;
		}

		s[i].lutime = 0.;
		s[i].r += gtime * s[i].gr;
	}

	r += gtime * growthrate;       // r defined at gtime = 0
	rtime += gtime;
	gtime = 0.;

	if (rescale == true)
		process(N);
}


//==============================================================
// Run time
//==============================================================
void Box::runTime() {
	time(&end);
	std::cout << "run time = " << difftime(end, start) << std::endl;
}


//Write configuration to output file
void Box::writeConfiguration(const char* wconfigfile, int iteration) {
	int count;

	if (gtime != 0.)   // Synchronize spheres if not currently synchronized
		synchronize(false);

	std::ofstream output;
	//if (iteration == 0)
		output = std::ofstream(wconfigfile); // Overwrite file if this is initial configuration
	//else
	//	output = std::ofstream(wconfigfile, std::ios_base::app); // Append new configuration to existing file

	//Formatting
	output.setf(std::ios::fixed, std::ios::floatfield);

	count = 0; // Number of spheres of first species
	for (int i = 0; i < N; i++) {
		if (s[i].species == 1) count += 1;
	}

	// make header
	output << "Iteration: " << iteration << std::endl;
	//output << DIM << std::endl;
	//output << N << " " << particle_sizes.size() << std::endl;
	//for (const auto& c : particle_count) {
	//	output << c << " ";
	//}
	//output << std::endl;

	//// Lattice vectors:
	//for (int i = 0; i < DIM; i++)
	//	for (int j = 0; j < DIM; j++)
	//	{
	//		if (i == j)
	//			output << 1 << " ";
	//		else
	//			output << 0 << " ";
	//	}
	//output << std::endl;

	//// Boundary conditions:
	//for (int i = 0; i < DIM; i++)
	//	output << "T" << " ";  // T = periodic BC
	//output << std::endl;

	//Write position and radius
	for (int i = 0; i < N; i++) {
		for (int k = 0; k < DIM; k++) {
			output << std::right << std::setw(20) << std::setprecision(16) << s[i].x[k] << " ";
		}
		output << std::right << std::setw(20) << std::setprecision(16) << s[i].r << std::endl;
	}

	output.close();
}
