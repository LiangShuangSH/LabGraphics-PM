#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include "RoomObject.h"
#include <fstream>
#include <iterator>
#include <random>

using namespace std;

int ACCESS_SPACE = 4;
double STANDARD_GAP = 0.2;
double NO_GAP = 0.0;
double SCALE = 0.08;

vector<RoomObject*> translate_furnitures(vector<RoomObject*> fnts, double temp, double walls[], int &chosen_index, int static_fnts);
vector<RoomObject*> rotate_furnitures(vector<RoomObject*> fnts, double temp, int &chosen_index, int static_fnts);
double get_accessibility(RoomObject* f1, RoomObject* f2, int k);
double cost_function(vector<RoomObject*> fnts, double walls[]);
double get_prior_distance_cost(FurnitureCode type, double distance);
double get_prior_theta_cost(FurnitureCode type, double theta);
void find_closest_wall(int& index, double& min_distance, double distance[], int size);
double get_theta(int index, double orientation);
FurnitureCode Hashit(string type);
void deep_copy_vector(vector<RoomObject*> v, vector<RoomObject*> &new_v);
double constrainAngle(double x);
double angleDiff(double a, double b);

int main(int argc, char** argv) {
	
	/*=================================*/
	/* Test Input */
	vector<RoomObject*> furnitures;

	/* Add floor and walls for testing */
	// Floor and Walls;
	RoomObject floor("Floor");
	floor.pos[1] = 0.01;
	furnitures.push_back(&floor);

	RoomObject wall_east("Wall1");
	wall_east.pos[0] = 3.0;
	wall_east.orientation = 0.0;
	wall_east.width = 6.0 * SCALE;
	wall_east.depth = 0.1 * SCALE;
	RoomObject wall_south("Wall1");
	wall_south.pos[2] = 3.0;
	wall_south.orientation = 1.5 * M_PI;
	wall_south.width = 6.0 * SCALE;
	wall_south.depth = 0.1 * SCALE;
	RoomObject wall_west("Wall1");
	wall_west.pos[0] = -3.0;
	wall_west.orientation = M_PI;
	wall_west.width = 6.0 * SCALE;
	wall_west.depth = 0.1 * SCALE;
	RoomObject wall_north("Wall2");
	wall_north.pos[2] = -3.0;
	wall_north.orientation = M_PI_2;
	wall_north.width = 6.0 * SCALE;
	wall_north.depth = 0.1 * SCALE;

	wall_east.update_access_space(NO_GAP);
	wall_south.update_access_space(NO_GAP);
	wall_west.update_access_space(NO_GAP);
	wall_north.update_access_space(NO_GAP);

	furnitures.push_back(&wall_east);
	furnitures.push_back(&wall_south);
	furnitures.push_back(&wall_west);
	furnitures.push_back(&wall_north);

	int STATIC_FNTS = 5;

	RoomObject sofa1("Sofa");
	sofa1.pos[0] = 0.0;
	sofa1.pos[2] = 0.0;
	sofa1.height = 0.3 * SCALE;
	sofa1.width = 6.0 * SCALE;
	sofa1.depth = 2.5 * SCALE;

	RoomObject sofa2("Sofa");
	sofa2.pos[0] = 0.0;
	sofa2.pos[2] = 0.0;
	sofa2.height = 0.3 * SCALE;
	sofa2.width = 6.0 * SCALE;
	sofa2.depth = 2.5 * SCALE;

	sofa1.update_access_space(STANDARD_GAP);
	sofa2.update_access_space(STANDARD_GAP);

	RoomObject sofa3("Sofa");
	sofa3.width = 6.0 * SCALE;
	sofa3.depth = 2.5 * SCALE;
	RoomObject sofa4("Bookcase");
	sofa4.width = 1.0 * SCALE;
	sofa4.depth = 0.4 * SCALE;
	RoomObject sofa5("Bookcase");
	sofa5.width = 1.0 * SCALE;
	sofa5.depth = 0.4 * SCALE;

	sofa3.update_access_space(STANDARD_GAP);
	sofa4.update_access_space(STANDARD_GAP);
	sofa5.update_access_space(STANDARD_GAP);

	furnitures.push_back(&sofa1);
	furnitures.push_back(&sofa2);
	furnitures.push_back(&sofa3);
	furnitures.push_back(&sofa4);
	furnitures.push_back(&sofa5);

	double walls[4];
	walls[0] = 3.0;
	walls[1] = 3.0;
	walls[2] = -3.0;
	walls[3] = -3.0;

	/*=================================*/

	//test cost
	//double test_cost = cost_function(furnitures, walls);

	// Simulated Annealing
	double curr_cost = cost_function(furnitures, walls);
	double new_cost = 0.0;;
	vector<RoomObject*> new_fnts;
	// Anealing Parameters
	int MAX_ITER = 10000000;
	double T0 = 10000.0;
	double T = T0; // Initial Temperature
	double beta = 0.01; // Boltzmann-like constant
	// Annealing starts
	for (int i = 0; i < MAX_ITER; i++) {
		// Pick a new layout for furnitures
		int chosen_index = 0;
		// First translate
		new_fnts = translate_furnitures(furnitures, T, walls, chosen_index, STATIC_FNTS);
		new_cost = cost_function(new_fnts, walls);
		double accept_prob = min(1.0, exp((1.0 /(beta * T)) * (curr_cost - new_cost)));
		srand(1);
		// Update position
		if (accept_prob >= ((double)rand() / (RAND_MAX))) {
			// furnitures = new_fnts;
			furnitures.at(chosen_index)->pos[0] = new_fnts.at(chosen_index)->pos[0];
			furnitures.at(chosen_index)->pos[2] = new_fnts.at(chosen_index)->pos[2];
			curr_cost = new_cost;
		}
		// Temperature drops
		T = T0 * pow(0.99, i);
		// Stop when T = 0.0
		if (T <= 0.01) {
			break;
		}
		cout << curr_cost << endl;

		// Then Rotate
		new_fnts = rotate_furnitures(furnitures, T, chosen_index, STATIC_FNTS);
		new_cost = cost_function(new_fnts, walls);
		accept_prob = min(1.0, exp((1.0 / (beta * T)) * (curr_cost - new_cost)));
		// Update Rotation
		if (accept_prob >= ((double)rand() / (RAND_MAX))) {
			furnitures.at(chosen_index)->orientation = new_fnts.at(chosen_index)->orientation;
			curr_cost = new_cost;
		}
		// Temperature drops
		T = T0 * pow(0.99, i);
		// Stop when T = 0.0
		if (T <= 0.01) {
			break;
		}
		cout << curr_cost << endl;
	}
	

	/* ============================= */

	// Output the furnitures
	ofstream out;
	out.open("furnitures.txt");
	for (int i = 0; i < furnitures.size(); i++) {
		RoomObject* ro = furnitures.at(i);
		out << "type=" << ro->type << "\n";
		out << "x=" << ro->pos[0] << "\n";
		out << "y=" << ro->pos[1] << "\n";
		out << "z=" << ro->pos[2] << "\n";
		out << "orientation=" << ro->orientation << "\n";
		out << "\n";
	}
	out.close();
	
	return 0;
}

vector<RoomObject*> translate_furnitures(vector<RoomObject*> fnts, double temp, double walls[], int &chosen_index, int static_fnts) {
	vector<RoomObject*> new_fnts;
	deep_copy_vector(fnts, new_fnts);
	double deviation = 0.05 * temp;
	
	random_device rd1;
	mt19937 generator(rd1());
	normal_distribution<double> distribution(0.0, deviation);
	double dx = distribution(generator);
	double dz = distribution(generator);

	random_device rd2; // obtain a random number from hardware
	mt19937 eng(rd2()); // seed the generator
	uniform_int_distribution<> distr(0 + static_fnts, new_fnts.size() - 1); // define the range

	chosen_index = distr(eng);
	RoomObject* the_chosen_one = new_fnts.at(chosen_index);
	the_chosen_one->pos[0] += dx;
	the_chosen_one->pos[0] = min(the_chosen_one->pos[0], walls[0]);
	the_chosen_one->pos[0] = max(the_chosen_one->pos[0], walls[2]);
	the_chosen_one->pos[2] += dz;
	the_chosen_one->pos[2] = min(the_chosen_one->pos[2], walls[1]);
	the_chosen_one->pos[2] = max(the_chosen_one->pos[2], walls[3]);
	the_chosen_one->update_access_space(STANDARD_GAP);

	return new_fnts;
}

vector<RoomObject*> rotate_furnitures(vector<RoomObject*> fnts, double temp, int &chosen_index, int static_fnts) {
	vector<RoomObject*> new_fnts;
	deep_copy_vector(fnts, new_fnts);
	double deviation = 0.01 * temp;

	random_device rd1;
	mt19937 generator(rd1());
	normal_distribution<double> distribution(0.0, deviation);
	double da = distribution(generator);

	random_device rd2; // obtain a random number from hardware
	mt19937 eng(rd2()); // seed the generator
	uniform_int_distribution<> distr(0 + static_fnts, new_fnts.size() - 1); // define the range

	chosen_index = distr(eng);
	RoomObject* the_chosen_one = new_fnts.at(chosen_index);
	the_chosen_one->orientation += da;
	the_chosen_one->orientation = fmod(the_chosen_one->orientation, 2.0 * M_PI);
	the_chosen_one->update_access_space(STANDARD_GAP);

	return new_fnts;
}

// fnts: current funitures layout
// new_fnts: new funitures layout
// walls[]: 4 walls, assume we are in a standard 4-wall house
double cost_function(vector<RoomObject*> fnts, double walls[]) {
	double cost = 0.0;
	double w_a = 0.1;
	double w_pd = 1.0;
	double w_td = 10.0;

	RoomObject* fnt1 = NULL;
	RoomObject* fnt2 = NULL;

	for (int i = 0; i < fnts.size(); i++) {
		fnt1 = fnts.at(i);
		// Calculate distance to walls
		double d[4] = {abs(walls[0] - fnt1->pos[0]),	// distance to east wall
						abs(walls[1] - fnt1->pos[2]),	// distance to south wall
						abs(walls[2] - fnt1->pos[0]),   // distance to west wall
						abs(walls[3] - fnt1->pos[2]) };	// distance to north wall

		double min_d, theta;
		int wall_index;
		find_closest_wall(wall_index, min_d, d, 4);
		theta = get_theta(wall_index, fnt1->orientation);

		// Prior knowledge about the furniture
		cost += w_pd * get_prior_distance_cost(Hashit(fnt1->type), min_d) + w_td * get_prior_theta_cost(Hashit(fnt1->type), theta);

		for (int j = 0; j < fnts.size(); j++) {
			fnt2 = fnts.at(j);
			for (int k = 0; k < ACCESS_SPACE; k++) {
				// Accessibility
				cost += w_a * get_accessibility(fnt1, fnt2, k);
			}
		}
	}

	return cost;
}

double get_accessibility(RoomObject* f1, RoomObject* f2, int k) {
	double result = 0.0;
	result = pow(f1->pos[0] - f2->access_space_center[k][0], 2) + pow(f1->pos[2] - f2->access_space_center[k][1], 2);
	result = sqrt(result);
	result = result / (f1->get_diagnol() + f2->access_space_diag[k]);
	result = 1 - result;
	result = max(0.0, result);
	return result;
}

void find_closest_wall(int& index, double& min_distance, double distance[], int size) {
	min_distance = distance[0];
	index = 0;
	for (int i = 1; i < size; i++) {
		if (min_distance > distance[i]) {
			min_distance = distance[i];
			index = i;
		}
	}
}

double get_theta(int index, double orientation) {
	double theta = 0.0;

	switch (index) {
	case 0: theta = angleDiff(orientation, 0.0);
			break;
	case 1: theta = angleDiff(orientation, (-M_PI_2));
			break;
	case 2: theta = angleDiff(orientation, M_PI);
			break;
	case 3: theta = angleDiff(orientation, M_PI_2);
			break;
	default: theta = 0.0;
			 break;
	}

	return theta;
}

FurnitureCode Hashit(string type) {
	if (type == "House") {
		return House;
	}
	else if (type == "Wall") {
		return Wall;
	}
	else if (type == "Floor") {
		return Floor;
	}
	else if (type == "Bed") {
		return Bed;
	}
	else if (type == "Sofa") {
		return Sofa;
	}
	else if (type == "Desk") {
		return Desk;
	}
	else if (type == "Chair") {
		return Chair;
	}
	else if (type == "Bookcase") {
		return Bookcase;
	}
	return Undefined;
}

void deep_copy_vector(vector<RoomObject*> v, vector<RoomObject*> &new_v) {
	for (int i = 0; i < v.size(); i++) {
		new_v.push_back(v.at(i)->clone());
	}
}

double constrainAngle(double x) {
	x = fmod(x + M_PI, 2.0*M_PI);
	if (x < 0.0)
		x += 2.0*M_PI;
	return x - M_PI;
}

double angleDiff(double a, double b) {
	double dif = fmod(b - a + M_PI, 2*M_PI);
	if (dif < 0.0)
		dif += 2*M_PI;
	return dif - M_PI;
}

/* ================================================================ */
/* Prior Knowledge */

double get_prior_distance_cost(FurnitureCode type, double distance) {
	double cost = 0.0;
	switch (type) {
	case Bed: cost = abs(distance - 0.5);
		break;
	case Sofa: cost = abs(distance - 0.5);
		break;
	case Desk: cost = abs(distance - 0.5);
		break;
	case Bookcase: cost = abs(distance - 0.5);
		break;
	default: cost = 0.0;
		break;
	}
	return cost;
}

double get_prior_theta_cost(FurnitureCode type, double theta) {
	double cost = 0.0;
	switch (type) {
	case Bed: cost = abs(angleDiff(theta, M_PI));
		break;
	case Sofa: cost = abs(angleDiff(theta, M_PI));
		break;
	case Desk: cost = abs(angleDiff(theta, M_PI));
		break;
	case Bookcase: cost = abs(angleDiff(theta, M_PI));
		break;
	default: cost = 0.0;
		break;
	}
	return cost;
}

/* ================================================================ */