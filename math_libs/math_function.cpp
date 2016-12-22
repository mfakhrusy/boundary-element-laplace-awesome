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

//green function -> for g calc
double Math_Function::green_f(double x, double x_ref, double y, double y_ref) {
	
	if (x == x_ref && y == y_ref) {
		return 0;
	} else {
		return log(pow(x - x_ref, 2) + pow(y - y_ref, 2));
	}
}

//derivation of green function in x (dG/dx)
double Math_Function::dgreen_f_dx(double x, double x_ref, double y, double y_ref) {
	
	if (x == x_ref) {
		return 0;
	} else {
		return 2*(x - x_ref)/(pow(x - x_ref,2) + pow(y - y_ref,2));
	}
}

//derivation of green function in y (dG/dy)
double Math_Function::dgreen_f_dy(double x, double x_ref, double y, double y_ref) {
	
	if (y == y_ref) {
		return 0;
	} else {
		return 2*(y - y_ref)/(pow(x - x_ref,2) + pow(y - y_ref,2));
	}
}

//derivation in n -> for h calc
double Math_Function::dgreen_f_dn(double x, double x_ref, double nx, double y, double y_ref, double ny) {
	
	return dgreen_f_dx(x, x_ref, y, y_ref)*nx + dgreen_f_dy(x, x_ref, y, y_ref)*ny;
}

double Math_Function::knonecker_delta(int i, int j) {
	
	if (i == j) {
		return 1;
	} else {
		return 0;
	}
}

double Math_Function::discriminant(double a, double b, double c) {
	
	return b*b - 4*a*c;
}

double Math_Function::integral_simpson(std::vector<double> x, std::vector<double> f_x) {
	
	int max_node	=	x.size();
	double x_i	=	x[0];
	double x_f	=	x[max_node - 1];
	double temp	=	(x_f - x_i)/(3.0*static_cast<double>(max_node));

	double sum_even = 0;
	for (auto i = 2; i < max_node-1; i = i+2) {
		sum_even = sum_even + f_x[i];
	}

	double sum_odd = 0;
	for (auto i = 1; i < max_node-1; i = i+2) {
		sum_odd = sum_odd + f_x[i];
	}

	double result	=	temp*(f_x[0] + f_x[max_node-1] + 4*sum_odd + 2*sum_even);
//	std::cout << "tututu: " << result << " " << temp << " " << f_x[0] << " " << f_x[max_node-1] << std::endl;
	return result;

}
