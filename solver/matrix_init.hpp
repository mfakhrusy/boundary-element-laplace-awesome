#ifndef MATRIX_INIT_HPP
#define MATRIX_INIT_HPP
#include "../math_libs/math_function.hpp"

class Math_Function;
class Matrix_Init {

	Math_Function math_f;

	double matrix_init_phi_calc(double x, double y, Parameters pars);
	double matrix_init_dphi_dn_calc(double nx, double ny, Parameters pars);
	std::vector<double> matrix_init_g1_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node);
	std::vector<double> matrix_init_g2_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node);
	std::vector<double> matrix_init_h1_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node);
	std::vector<double> matrix_init_h2_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node);
	std::vector<double> matrix_init_H_rhs_calc(Airfoil_Parameters airfoil_pars, Parameters pars);
	std::vector<double> matrix_init_H_lhs_calc(Airfoil_Parameters airfoil_pars, Parameters pars);

	std::vector<double> matrix_init_rhs_matrix_calc(Airfoil_Parameters airfoil_pars, Parameters pars, Variables &vars);
	std::vector<std::vector<double>> matrix_init_lhs_matrix_calc(Airfoil_Parameters airfoil_pars, Parameters pars, Variables &vars);

	public:
		void matrix_init_main_computation(Airfoil_Parameters airfoil_pars, Parameters pars, Variables &vars);
};

#endif
