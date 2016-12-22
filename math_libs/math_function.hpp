#ifndef MATH_FUNCTION_HPP
#define MATH_FUNCTION_HPP

class Math_Function {

	double length_two_points(double x1, double x2, double y1, double y2);
	double angle_cos_rule(double A_this, double B, double C);
	double avg_n(double n_1, double n_2, double n_couple_1, double n_couple_2); 
	double green_f(double x, double x_ref, double y, double y_ref);
	double dgreen_f_dx(double x, double x_ref, double y, double y_ref);
	double dgreen_f_dy(double x, double x_ref, double y, double y_ref);
	double dgreen_f_dn(double x, double x_ref, double nx, double y, double y_ref, double ny);
	double knonecker_delta(int i, int j);
	double discriminant(double a, double b, double c);
	double integral_simpson(std::vector<double> x, std::vector<double> f_x);

	//friend classes below can access the functions above (which are privates)
	friend class Airfoil;
	friend class Matrix_Init;
	friend class Matrix_Solver;
};

#endif
