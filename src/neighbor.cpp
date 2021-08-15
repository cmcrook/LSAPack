#include "vec.hpp"
#include <iostream>
#include "box.hpp"

//  Class neighbor
Neighbor::Neighbor(int i_i) : i(i_i) {}

//  Class collision
Collision::Collision(int i_i, Box* b_i)
	: Neighbor(i_i),
	box(b_i) {
	ctime = dblINF;
	cpartner = i;
}

// Operation is finding the next collision from a given cell
void Collision::Operation(int j, vector<DIM, int>& pboffset) {
	box->predictCollision(i, j, pboffset, ctime, cpartner, cpartnerpboffset);
}

