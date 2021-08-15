#pragma once

#include "vec.hpp"

class Sphere {
 public:

  // constructor and destructor
  Sphere();
  Sphere(const Sphere& s);
  Sphere(int i_i, vector<DIM> x, vector<DIM, int> cell_i, double lutime_i, 
	 double r_i, double gr_i, double m_i, int species_i);
  ~Sphere();

 //variables
  int id;                          // sphere ID

  // impending event
  Event nextevent;                // next event...can be collision or transfer
  Event nextcollision;            // next collision if next event is transfer
  // maybe nextnext event
  
  // past information
  double lutime;                  // last update time
  vector<DIM, int> cell;          // cell that it belongs to
  vector<DIM, double> x;          // position
  vector<DIM, double> v;          // velocity
  double r;                       // sphere radius
  double gr;                      // sphere growth rate
  double m;                       // sphere mass
  int species;                    // species number (not used during the MD)
  // make sure efficent in memory
};