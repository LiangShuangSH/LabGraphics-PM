#pragma once

#include <iostream>
#include <vector>

using namespace std;

enum FurnitureCode {
	Undefined = -1,
	House = 0,
	Wall = 1,
	Floor = 2,
	Bed = 3,
	Sofa = 4,
	Desk = 5,
	Chair = 6,
	Bookcase = 7
};

class RoomObject {
public:
	string type; // such as "bed", "chair"
	double pos[3] = { 0.0 }; // contains 3 float: x. y, z
	double orientation = 0.0; // Here we only consider the yaw of an object
	double height = 0.0, width = 0.0, depth = 0.0; // The size of object
	double access_space_center[4][2];
	double access_space_diag[4];

	RoomObject(string type) { this->type = type; }
	double get_diagnol();
	void update_access_space(double gap);
	RoomObject* clone();

private:
	RoomObject* _prev = NULL;
	vector<RoomObject*> _next;
};