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
#include <ctime>

using namespace std;

int ACCESS_SPACE = 4;
double STANDARD_GAP = 0.0;
double NO_GAP = 0.0;
int STATIC_FNTS = 5;

vector<RoomObject*> translate_furnitures(vector<RoomObject*> fnts, double temp, double walls[], int &chosen_index, int static_fnts);
vector<RoomObject*> rotate_furnitures(vector<RoomObject*> fnts, double temp, int &chosen_index, int static_fnts);
double get_accessibility(RoomObject* f1, RoomObject* f2, int k);
double cost_function(vector<RoomObject*> fnts, double walls[]);
double get_prior_distance_cost(FurnitureCode type, double distance);
double get_prior_theta_cost(FurnitureCode type, double theta);
void find_closest_wall(int& index, double& min_distance, double distance[], int size);
double get_theta(int index, double orientation);
double get_out_of_wall_punishment(RoomObject* f, double walls[]);
FurnitureCode Hashit(string type);
void deep_copy_vector(vector<RoomObject*> v, vector<RoomObject*> &new_v);
double constrainAngle(double x);
double angleDiff(double a, double b);
void visualize_in_matlab_code(vector<RoomObject*> fnts, string matlab_file);
vector<double*> find_four_corners(RoomObject* f);
double get_pairwise_distance_cost(RoomObject* f);
double get_pairwise_theta_cost(RoomObject* f);
double get_pairwise_side_cost(RoomObject* f);
int find_index(RoomObject* f, vector<RoomObject*> fnts);
void random_start(vector<RoomObject*> &fnts, double walls[]);
void record_cost(vector<double> costs, string out_file);

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
	wall_east.pos[0] = 4.0;
	wall_east.orientation = 0.0;
	wall_east.width = 8.0;
	wall_east.depth = 0.1;
	RoomObject wall_south("Wall1");
	wall_south.pos[2] = 4.0;
	wall_south.orientation = -M_PI_2;
	wall_south.width = 8.0;
	wall_south.depth = 0.1;
	RoomObject wall_west("Wall1");
	wall_west.pos[0] = -4.0;
	wall_west.orientation = M_PI;
	wall_west.width = 8.0;
	wall_west.depth = 0.1;
	RoomObject wall_north("Wall2");
	wall_north.pos[2] = -4.0;
	wall_north.orientation = M_PI_2;
	wall_north.width = 8.0;
	wall_north.depth = 0.1;

	wall_east.update_access_space(NO_GAP);
	wall_south.update_access_space(NO_GAP);
	wall_west.update_access_space(NO_GAP);
	wall_north.update_access_space(NO_GAP);

	furnitures.push_back(&wall_east);
	furnitures.push_back(&wall_south);
	furnitures.push_back(&wall_west);
	furnitures.push_back(&wall_north);

	RoomObject sofa1("Sofa");
	sofa1.pos[0] = 0.0;
	sofa1.pos[2] = 0.0;
	sofa1.height = 0.3;
	sofa1.width = 2.0;
	sofa1.depth = 1.0;

	RoomObject sofa2("Bed");
	sofa2.pos[0] = 0.0;
	sofa2.pos[2] = 0.0;
	sofa2.height = 0.7;
	sofa2.width = 1.0;
	sofa2.depth = 1.5;

	sofa1.update_access_space(STANDARD_GAP);
	sofa2.update_access_space(STANDARD_GAP);

	RoomObject sofa3("PCTable");
	sofa3.width = 1.5;
	sofa3.depth = 0.6;
	RoomObject sofa4("Chair");
	sofa4.width = 0.5;
	sofa4.depth = 0.6;
	sofa4.pairwise = "Chair and PC Table";
	RoomObject sofa5("Bookcase");
	sofa5.width = 1.0;
	sofa5.depth = 0.4;
	RoomObject f6("Sofa");
	f6.pos[0] = 0.0;
	f6.pos[2] = 0.0;
	f6.height = 0.3;
	f6.width = 2.0;
	f6.depth = 1.0;
	

	sofa3.update_access_space(STANDARD_GAP);
	sofa4.update_access_space(STANDARD_GAP);
	sofa5.update_access_space(STANDARD_GAP);
	f6.update_access_space(STANDARD_GAP);

	sofa1.set_prev_node(&floor);
	sofa2.set_prev_node(&floor);
	sofa3.set_prev_node(&floor);
	sofa4.set_prev_node(&sofa3);
	sofa5.set_prev_node(&floor);
	f6.set_prev_node(&floor);

	furnitures.push_back(&sofa1);
	furnitures.push_back(&sofa2);
	furnitures.push_back(&sofa3);
	furnitures.push_back(&sofa4);
	furnitures.push_back(&sofa5);
	furnitures.push_back(&f6);

	double walls[4];
	walls[0] = 4.0;
	walls[1] = 4.0;
	walls[2] = -4.0;
	walls[3] = -4.0;

	/*=================================*/
	//test cost
	//double test_cost = cost_function(furnitures, walls);

	// Simulated Annealing
	random_start(furnitures, walls);
	double curr_cost = cost_function(furnitures, walls);
	double new_cost = 0.0;
	vector<RoomObject*> new_fnts;
	// Anealing Parameters
	int MAX_ITER = 10000000;
	double T0 = 10000.0;
	double T = T0; // Initial Temperature
	double beta = 0.03; // Boltzmann-like constant
	vector<double> costs;
	// Annealing starts
	for (int i = 0; i < MAX_ITER; i++) {
		// Pick a new layout for furnitures
		int chosen_index = 0;
		// First translate
		new_fnts = translate_furnitures(furnitures, T, walls, chosen_index, STATIC_FNTS);
		new_cost = cost_function(new_fnts, walls);
		double accept_prob = min(1.0, exp((1.0 /(beta * T)) * (curr_cost - new_cost)));

		// Update position
		double rand_num = (double)rand() / RAND_MAX;
		if (accept_prob >= rand_num) {
			furnitures = new_fnts;
			curr_cost = new_cost;
		}
		costs.push_back(curr_cost);
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
			furnitures = new_fnts;
			//furnitures.at(chosen_index)->orientation = new_fnts.at(chosen_index)->orientation;
			curr_cost = new_cost;
		}
		costs.push_back(curr_cost);
		// Temperature drops
		T = T0 * pow(0.99, i);
		// Stop when T = 0.0
		if (T <= 0.01) {
			break;
		}
		cout << curr_cost << endl;

		if ((i % 100) == 0) {
			string file = "Matlab" + to_string(i) + ".m";
			visualize_in_matlab_code(furnitures, file);
		}
	}
	record_cost(costs, "costs.txt");

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

	visualize_in_matlab_code(furnitures, "Mat_Visual.m");
	
	return 0;
}

vector<RoomObject*> translate_furnitures(vector<RoomObject*> fnts, double temp, double walls[], int &chosen_index, int static_fnts) {
	vector<RoomObject*> new_fnts;
	deep_copy_vector(fnts, new_fnts);
	double deviation = max(0.001 * temp, 0.5);
	
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
	double w_a = 7.0;
	double w_pd = 5.0;
	double w_td = 10.0;
	double w_wall_punish = 1.0;
	double w_pr_dis = 1.0;
	double w_pr_theta = 10.0;
	double w_pr_side = 5.0;

	RoomObject* fnt1 = NULL;
	RoomObject* fnt2 = NULL;

	for (int i = STATIC_FNTS; i < fnts.size(); i++) {
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
		cost += w_wall_punish * get_out_of_wall_punishment(fnt1, walls);
		cost += w_pr_dis * get_pairwise_distance_cost(fnt1);
		cost += w_pr_theta * get_pairwise_theta_cost(fnt1);
		cost += w_pr_side * get_pairwise_side_cost(fnt1);

		for (int j = STATIC_FNTS; j < fnts.size(); j++) {
			fnt2 = fnts.at(j);
			for (int k = 0; k < ACCESS_SPACE; k++) {
				// Accessibility
				cost += w_a * get_accessibility(fnt1, fnt2, k);
			}
		}
	}

	return cost;
}

void random_start(vector<RoomObject*> &fnts, double walls[]) {
	double right = walls[0];
	double bottom = walls[1];
	double left = walls[2];
	double up = walls[3];

	RoomObject* f;
	double random;
	srand(time(NULL));

	for (int i = STATIC_FNTS; i < fnts.size(); i++) {
		f = fnts.at(i);	

		random = (double)rand() / RAND_MAX;
		f->pos[0] = left + random * (right - left);

		random = (double)rand() / RAND_MAX;
		f->pos[2] = up + random * (bottom - up);

		random = (double)rand() / RAND_MAX;
		f->orientation = -M_PI + random * 2.0 * M_PI;
	}
}

double get_accessibility(RoomObject* f1, RoomObject* f2, int k) {
	double result = 0.0;

	result = pow(f1->pos[0] - f2->access_space_center[k][0], 2) + pow(f1->pos[2] - f2->access_space_center[k][1], 2);
	result = sqrt(result);
	result = result / (0.5 * (f1->get_diagnol() + f2->access_space_diag[k]));
	if (result > 1.0) {
		return 0.0;
	}
	result = 1.0 - pow(result, 0.5);
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

double get_out_of_wall_punishment(RoomObject* f, double walls[]) {
	double punishment = 100.0;
	vector<double*> vertice = find_four_corners(f);
	double x, y;
	for (int i = 0; i < vertice.size(); i++) {
		x = vertice.at(i)[0];
		y = vertice.at(i)[1];

		if ((x > walls[0]) || (x < walls[2])) {
			return punishment;
		}
		else if ((y > walls[1]) || (y < walls[3])) {
			return punishment;
		}
	}
	return 0.0;
}

void visualize_in_matlab_code(vector<RoomObject*> fnts, string matlab_file) {
	ofstream output;
	output.open(matlab_file);

	string code = "";
	RoomObject* curr;

	code += "clear;\n";
	code += "figure;\n";
	code += "hold on;\n";
	for (int i = 1; i < fnts.size(); i++) {
		curr = fnts.at(i);
		string x1, x2, y1, y2, rotation, c_x, c_z, color;
		string idx = to_string(i);
		x1 = to_string(-(curr->width / 2.0));
		x2 = to_string(curr->width / 2.0);
		y1 = to_string(-(curr->depth / 2.0));
		y2 = to_string(curr->depth / 2.0);
		rotation = to_string(curr->orientation + M_PI_2);
		c_x = to_string(curr->pos[0]);
		c_z = to_string(curr->pos[2]);
		if (curr->type.find("Wall") != string::npos) {
			color = "\'r-\'";
		}
		else if (curr->type.find("PCTable") != string::npos) {
			color = "\'g-\'";
		}
		else if (curr->type.find("Chair") != string::npos) {
			color = "\'g-\'";
		}
		else {
			color = "\'b-\'"; // blue
		}

		code += "X" + idx + " = [" + x1 + "," + x2 + "," + x2 + "," + x1 + "," + x1 + "];\n"; // Xi = [x1, x2, x2, x1, x1];
		code += "Y" + idx + " = [" + y1 + "," + y1 + "," + y2 + "," + y2 + "," + y1 + "];\n";   // Yi = [y1, y1, y2, y2, y1];
		code += "h" + idx + " = plot(X" + to_string(i) + ", Y" + to_string(i) + ", " + color + ", 'LineWidth', 3);\n";
		code += "t" + idx + " = hgtransform('Parent', gca);\n"; // t = hgtransform('Parent',gca);
		code += "set(h" + idx + ", \'Parent\', t" + idx + ")\n"; // set(h, 'Parent', t)
		code += "Txy" + idx + " = makehgtform(\'zrotate\', " + rotation + ");\n"; // Txy = makehgtform('zrotate', rotation);
		code += "Rxy" + idx + " = makehgtform(\'translate\', [" + c_x + " " + c_z + " 0]);\n";// Rxy = makehgtform('translate', [pos[0] pos[2] 0]);
		code += "Trans" + idx + " = " + "Rxy" + idx + " * " + "Txy" + idx + ";\n"; // Trans = Rxy * Txy;
		code += "set(t" + idx + ",\'Matrix\', Trans" + idx + ")\n"; // set(t,'Matrix',Trans)
		code += "\n";
	}
	code += "set(gca,'Ydir','reverse')\n";
	code.append("xlim([-4, 4]);\n");
	code.append("ylim([-4, 4]);\n");
	code.append("axis square;\n");

	output << code;
	output.close();
}

void record_cost(vector<double> costs, string out_file) {
	ofstream out;
	out.open(out_file);
	for (int i = 0; i < costs.size(); i++) {
		out << costs.at(i) << "\n";
	}
	out.close();
}

vector<double*> find_four_corners(RoomObject* f) {
	vector<double*> corners;

	double sign_x, sign_y;
	double *coord;
	for (int i = 0; i < 4; i++) {
		if (i == 0 || i == 3) {
			sign_x = 1.0;
		}
		else {
			sign_x = -1.0;
		}

		if (i == 0 || i == 1) {
			sign_y = 1.0;
		}
		else {
			sign_y = -1.0;
		}
		coord = new double[2];
		coord[0] = f->pos[0] + (sign_x * f->depth / 2.0) * cos(f->orientation) - (sign_y * f->width / 2.0) * sin(f->orientation);
		coord[1] = f->pos[2] + (sign_x * f->depth / 2.0) * sin(f->orientation) + (sign_y * f->width / 2.0) * cos(f->orientation);
		corners.push_back(coord);
	}
	return corners;
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
	case PCTable: cost = abs(distance - 0.3);
		break;
	case Bookcase: cost = abs(distance - 0.25);
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
	case PCTable: cost = abs(angleDiff(theta, M_PI));
		break;
	case Bookcase: cost = abs(angleDiff(theta, M_PI));
		break;
	default: cost = 0.0;
		break;
	}
	return cost;
}

double get_pairwise_distance_cost(RoomObject* f) {
	double cost = 0.0;
	RoomObject* partner = f->get_prev_node();

	double distance = pow(partner->pos[0] - f->pos[0], 2) + pow(partner->pos[2] - f->pos[2], 2);
	distance = sqrt(distance);

	if (f->pairwise == "Chair and PC Table") {
		cost = abs(distance - 1.0);
	}
	else {
		cost = 0.0;
	}
	return cost;
}

// Defines the orientation relationship between partners
// For example the chair should face to desk
double get_pairwise_theta_cost(RoomObject* f) {
	double cost = 0.0;
	RoomObject* partner = f->get_prev_node();

	double theta_diff = angleDiff(f->orientation, partner->orientation);

	if (f->pairwise == "Chair and PC Table") {
		cost = abs(angleDiff(theta_diff, M_PI));
	}
	else {
		cost = 0.0;
	}
	return cost;
}

// Defines which side to put the object, like front, back or sides
double get_pairwise_side_cost(RoomObject* f) {
	double cost = 0.0;
	RoomObject* partner = f->get_prev_node();

	double distance = pow(partner->pos[0] - f->pos[0], 2) + pow(partner->pos[2] - f->pos[2], 2);
	distance = sqrt(distance);
	if (distance == 0.0) {
		return M_PI;
	}

	double pos_angle = asin(abs(f->pos[2] - partner->pos[2]) / distance);
	pos_angle = constrainAngle(pos_angle);
	double angle_diff = angleDiff(pos_angle, partner->orientation);
	
	if (f->pairwise == "Chair and PC Table") {
		cost = abs(angleDiff(angle_diff, 0.0));
	}
	else {
		cost = 0.0;
	}
	return cost;
}

/* ================================================================ */

/* ================================================================ */
/* Utility */
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
	else if (type == "PCTable") {
		return PCTable;
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
	// Copy Object
	for (int i = 0; i < v.size(); i++) {
		new_v.push_back(v.at(i)->clone());
	}
	// Re-link objects
	RoomObject* prev = NULL;
	int idx = 0;
	RoomObject* new_mother;
	RoomObject* new_child;
	for (int i = 0; i < v.size(); i++) {
		prev = v.at(i)->get_prev_node();
		if (prev != NULL) {
			idx = find_index(prev, v);
			new_mother = new_v.at(idx);
			new_child = new_v.at(i);
			new_child->set_prev_node(new_mother);
			new_mother->add_child(new_child);
		}
	}
}

double constrainAngle(double x) {
	x = fmod(x + M_PI, 2.0*M_PI);
	if (x < 0.0)
		x += 2.0*M_PI;
	return x - M_PI;
}

double angleDiff(double a, double b) {
	double dif = fmod(b - a + M_PI, 2 * M_PI);
	if (dif < 0.0)
		dif += 2 * M_PI;
	return dif - M_PI;
}

int find_index(RoomObject* f, vector<RoomObject*> fnts) {
	int idx = 0;
	for (int i = 0; i < fnts.size(); i++) {
		if (f == fnts.at(i)) {
			idx = i;
		}
	}
	return idx;
}

/* ================================================================ */