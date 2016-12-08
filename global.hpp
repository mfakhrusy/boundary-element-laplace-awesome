#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

struct Parameters {

	double angle_of_attack;
	double v_infinity;
	int max_node;
};

struct Airfoil_Parameters {
	
	std::vector<double> x;		//node position on x 
	std::vector<double> y;		//node position on y
	std::vector<double> s;		//element length
	std::vector<double> beta;	//beta is the angle between 2 adjacent tangent line
};



#endif
