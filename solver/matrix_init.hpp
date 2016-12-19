#ifndef MATRIX_INIT_HPP
#define MATRIX_INIT_HPP
#include "../math_libs/math_function.hpp"
#include "../misc/miscellaneous.hpp"

class Math_Function;
class Miscellaneous;
class Matrix_Init {

	Math_Function math_f;
	Miscellaneous misc;

	double matrix_init_phi_calc(double x, double y, Parameters pars);
	double matrix_init_dphi_dn_calc(double nx, double ny, Parameters pars);
	std::vector<double> matrix_init_g1_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node);
	std::vector<double> matrix_init_g2_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node);
	std::vector<double> matrix_init_h1_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node);
	std::vector<double> matrix_init_h2_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node);
	std::vector<double> matrix_init_H_rhs_calc(Airfoil_Parameters airfoil_pars, Parameters pars);
	std::vector<double> matrix_init_H_lhs_calc(Airfoil_Parameters airfoil_pars, Parameters pars);

	std::vector<double> matrix_init_rhs_matrix_calc(Airfoil_Parameters airfoil_pars, Parameters pars, Variables vars);
	std::vector<std::vector<double>> matrix_init_lhs_matrix_calc(Airfoil_Parameters airfoil_pars, Parameters pars, Variables vars);


	//================= 2nd version of g1 and g2 calculation ==============
	class G_Calc {

		Math_Function math_f;
		double a_calc(double x_ref, double x_j, double x_j1, double y_ref, double y_j, double y_j1);
		double b_calc(double x_ref, double x_j, double x_j1, double y_ref, double y_j, double y_j1);
		double c_calc(double x_ref, double x_j, double x_j1, double y_ref, double y_j, double y_j1);

		double I_1_calc(double a, double b, double c);
		double I_2_calc(double a, double b, double c);
		double I_3_calc(double a, double b, double c, double I_1, double I_2);
		double I_4_calc(double a, double b, double c, double I_3);

		public:
			//double g_ij_calc(double s, double I_3, double I_4);
			double g_ij_calc(double s, double x_ref, double x_j, double x_j1, double y_ref, double y_j, double y_j1);

	}g_calc;

	//lhs v2
	std::vector<std::vector<double>> matrix_init_lhs_matrix_calc_v2(Airfoil_Parameters airfoil_pars, Parameters pars, Variables vars);

	public:
		void matrix_init_main_computation(Airfoil_Parameters airfoil_pars, Parameters pars, Variables &vars);
};

#endif
