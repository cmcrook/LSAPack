#include "event.hpp" 

//==============================================================
//==============================================================
//  Class Event 
//==============================================================
//==============================================================


//==============================================================
// Constructor
//==============================================================
Event::Event(double time_i, int i_i, int j_i, vector<DIM,int> v_i):
  time(time_i),
  i(i_i),
  j(j_i),
  v(v_i)
{
}

Event::Event(double time_i, int i_i, int j_i):
  time(time_i),
  i(i_i),
  j(j_i)
{
}

Event::Event(const Event& e)
{
  time = e.time;
  i = e.i;
  j = e.j;
  v = e.v;
}

Event::Event()
{
}

  
//==============================================================
// Destructor
//==============================================================
Event::~Event() 
{
}

void Event::erase()
{
  time = dblINF;
  i = 0;
  j = 0;
}

bool Event::operator<(const Event& e) const
{
  return e.time < time;
}

bool Event::operator>(const Event& e) const
{
  return e.time > time;
}

 
