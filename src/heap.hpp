#pragma once
//---------------------------------------------------------------------------
// Event heap maker
//---------------------------------------------------------------------------

#include "event.hpp"
#include "sphere.hpp"

#include <vector>

class Heap {

 public:

  // constructor and destructor
  Heap(int maxsize, std::vector<Sphere>& s);
  Heap(const Heap &h);

  // variables
  int maxsize;   // max allowed number of events
  int N;         // current number of events
  std::vector<int> a;
  std::vector<Sphere>& s;
  std::vector<int> index;     // array of indices for each sphere
  //event minevent;


  // functions which operate on a binary heap
  
  void upheap(int k);
  void downheap(int k);
  void insert(int i);
  void replace(int i);
  int search(int j);
  void change(int i); 
  int extractmax();
  void print();
  void checkindex();
  
};