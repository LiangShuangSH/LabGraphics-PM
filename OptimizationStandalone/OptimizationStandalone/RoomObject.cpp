#include "RoomObject.h"
#define _USE_MATH_DEFINES
#include <math.h>

double RoomObject::get_diagnol() {
	double diagnol = 0.0;
	diagnol = pow(this->width, 2) + pow(this->depth, 2);
	diagnol = sqrt(diagnol);
	return diagnol;
}

void RoomObject::update_access_space(double gap) {
	// Update diagnol
	this->access_space_diag[0] = sqrt(pow(this->width, 2) + pow(gap + 0.5*this->depth, 2));
	this->access_space_diag[2] = sqrt(pow(this->width, 2) + pow(gap + 0.5*this->depth, 2));
	this->access_space_diag[1] = sqrt(pow(this->depth, 2) + pow(gap + 0.5*this->width, 2));
	this->access_space_diag[3] = sqrt(pow(this->depth, 2) + pow(gap + 0.5*this->width, 2));

	// Update 4 centers
	// 0:front 1:right 2:back 3:left
	for (int i = 0; i < 4; i++) {
		if (i % 2 == 0) {
			this->access_space_center[i][0] = this->pos[0] + 0.5 * (0.5 * depth + gap) * cos(orientation + i * M_PI_2);
			this->access_space_center[i][1] = this->pos[2] + 0.5 * (0.5 * depth + gap) * sin(orientation + i * M_PI_2);
		}
		else {
			this->access_space_center[i][0] = this->pos[0] + 0.5 * (0.5 * width + gap) * cos(orientation + i * M_PI_2);
			this->access_space_center[i][1] = this->pos[2] + 0.5 * (0.5 * width + gap) * sin(orientation + i * M_PI_2);
		}
	}
}

RoomObject* RoomObject::clone() {
	RoomObject* theClone = new RoomObject(this->type);

	theClone->pos[0] = this->pos[0];
	theClone->pos[1] = this->pos[1];
	theClone->pos[2] = this->pos[2];

	theClone->orientation = this->orientation;

	theClone->height = this->height;
	theClone->width = this->width;
	theClone->depth = this->depth;

	return theClone;
}