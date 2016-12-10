#ifndef MATH_FUNCTION_HPP
#define MATH_FUNCTION_HPP

class Math_Function {

	double length_two_points(double x1, double x2, double y1, double y2);
	double angle_cos_rule(double A_this, double B, double C);
	double avg_n(double n_1, double n_2, double n_couple_1, double n_couple_2); 

	//friend classes below can access the functions above (which are privates)
	friend class Airfoil;
};

#endif
