#pragma once

#include <iostream>
#include "RoomObject.h"

using namespace std;

struct Pairwise {
	RoomObject* origin;
	vector<string> partners;
	string relation;
};