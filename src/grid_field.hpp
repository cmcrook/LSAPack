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
#include "vec.hpp"

#include <vector>

// ======================================================================
// grid_field 
// ======================================================================

// A field of V-vectors on a D dimensional manifold

template<int D, class T>
class GridField {
  
 public:
  int elements;

 private:
  std::vector<T> f;
  vector<D, int> size;           // number of grid points for each dimension
  vector<D, int> offset;
 
 public:

  GridField();
  GridField(const vector<D, int>&);

  T& get(const vector<D, int>&);

  vector<D, int> get_size() const;
  void set_size(const vector<D, int>&);
  void set_size(const int);
  void initialize(const int i);
};


// grid_field
// ~~~~~~~~~~~~
template<int D, class T>
GridField<D, T>::GridField()
  : elements(0)
{
}


// grid_field
// ~~~~~~~~~~~~
template<int D, class T>
GridField<D, T>::GridField(const vector<D, int>& s)
{
  set_size(s);
}

// get_size
// ~~~~~~~~
template<int D, class T>
inline vector<D, int> GridField<D, T>::get_size() const
{
  return size;
}


// set_size
// ~~~~~~~~
template<int D, class T>
void GridField<D, T>::set_size(const vector<D, int>& s)
{
   f.clear();

  size = s;

  elements = 1;
  for(int i=0; i<D; i++) {
    offset[i] = elements;
    elements *= size.x[i];
  }

  f.resize(elements);
}


// set_size
// ~~~~~~~~
template<int D, class T>
void GridField<D, T>::set_size(const int s)
{
  vector<D, int> square;

  for(int k=0; k<D; k++)
    square[k] = s;

  set_size(square);
}


// get
// ~~~
template<int D, class T>
inline T& GridField<D, T>::get(const vector<D, int>& pos)
{
  int p=0;
  for(int i=0; i<D; i++)
    p += pos.x[i]*offset[i];

  if (p >= elements) {
      std::cout << "Elements: " << elements << std::endl;
      std::cout << "Size: " << size.x[0] << " " << size.x[1] << " " << size.x[2] << std::endl;
      std::cout << "Pos: " << pos.x[0] << " " << pos.x[1] << " " << pos.x[2] << std::endl;
      std::cout << "Offsets: " << offset[0] << " " << offset[1] << " " << offset[2] << std::endl;
      std::cout << "Index: " << p << std::endl;
      exit(-1);
  }


  return f[p];
}


// initialize
// ~~~
template<int D, class T>
void GridField<D, T>::initialize(const int value)
{
  for(int i=0; i<elements; i++)
    f[i] = value;
}
