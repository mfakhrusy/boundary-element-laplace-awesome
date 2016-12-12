#ifndef MATRIX_INIT_HPP
#define MATRIX_INIT_HPP
#include "../math_libs/math_function.hpp"

class Math_Function;
class Matrix_Init {

	Math_Function math_f;

	double matrix_init_calc(double nx, double ny, Parameters pars);
	std::vector<double> matrix_init_g1_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node);
	std::vector<double> matrix_init_g2_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node);
	std::vector<double> matrix_init_h1_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node);
	std::vector<double> matrix_init_h2_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node);

	public:
		void matrix_init_main_computation(Airfoil_Parameters airfoil_pars, Parameters pars, Variables &vars);
};

#endif
