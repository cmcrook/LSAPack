#include "box.hpp"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "vec.hpp"

//==============================================================
//==============================================================
//  Class Sphere: 
//==============================================================
//==============================================================


//==============================================================
// Constructor
//==============================================================
Sphere::Sphere()
{
}


//==============================================================
// Constructor
//==============================================================
Sphere::Sphere(const Sphere& s)
{
  id = s.id;
  x = s.x;
  v = s.v;
  cell = s.cell;
  lutime = s.lutime;
  nextevent = s.nextevent;
  nextcollision = s.nextcollision;
  r = s.r;
  gr = s.gr;
  m = s.m;
  species=s.species;
}

//==============================================================
// Constructor
//==============================================================
Sphere::Sphere(int i_i, vector<DIM> x_i, vector<DIM, int> cell_i, 
	       double lutime_i, double r_i, double gr_i, double m_i, int species_i):
  id(i_i),
  x(x_i),
  cell(cell_i),
  lutime(lutime_i),
  r(r_i),
  gr(gr_i),
  m(m_i),
  species(species_i)
{
}

//==============================================================
// Destructor
//==============================================================
Sphere::~Sphere() 
{

}

 
