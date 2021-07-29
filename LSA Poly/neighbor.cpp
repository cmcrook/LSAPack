#include "vec.hpp"
#include <iostream>
#include "box.hpp"

//  Class neighbor
neighbor::neighbor(int i_i) : i(i_i) {}

//  Class collision
Collision::Collision(int i_i, Box* b_i)
	: neighbor(i_i),
	b(b_i) {
	ctime = dblINF;
	cpartner = i;
}

// Operation is finding the next collision from a given cell
void Collision::Operation(int j, vector<DIM, int>& pboffset) {
	b->predictCollision(i, j, pboffset, ctime, cpartner, cpartnerpboffset);
}

