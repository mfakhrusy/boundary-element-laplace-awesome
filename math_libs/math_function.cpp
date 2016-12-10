#include "../global.hpp"
#include "math_function.hpp"

double Math_Function::length_two_points(double x1, double x2, double y1, double y2) {

	double delta_x	=	x2 - x1;
	double delta_y	=	y2 - y1;
	double length	=	sqrt(pow(delta_x,2) + pow(delta_y,2));

	return length;
}

double Math_Function::angle_cos_rule(double A_this, double B, double C) {

	double temp_cos	=	(pow(B,2) + pow(C,2) - pow(A_this,2))/(2*B*C);

	return acos(temp_cos);
}

double Math_Function::avg_n(double n_1, double n_2, double n_couple_1, double n_couple_2) {
	
	return (n_1 + n_2)/(sqrt(pow(n_1 + n_2, 2) + pow(n_couple_1 + n_couple_2,2)));

}
