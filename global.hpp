#ifndef GLOBAL_HPP
#define GLOBAL_HPP

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <numeric>	//performing accumulation/summation of vectors

struct Parameters {

	double v_freestream;	// v_infinity -> input
	double angle_of_attack;
	int max_node;
};

struct Airfoil_Parameters {
	
	std::vector<double> x;		//node position on x 
	std::vector<double> y;		//node position on y
	std::vector<double> s;		//element length
	std::vector<double> beta;	//beta is the angle between 2 adjacent tangent line
	std::vector<double> nx;		//normal derivative - x direction
	std::vector<double> ny;		//normal derivative - y direction
	std::vector<double> eta;	//shape functions
	std::vector<double> N1;
	std::vector<double> N2;
};

struct Variables {
	
	std::vector<double> rhs_matrix;
	std::vector<std::vector<double>> lhs_matrix;
};


#endif
