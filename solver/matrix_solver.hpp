#ifndef MATRIX_SOLVER_HPP
#define MATRIX_SOLVER_HPP
#include "../math_libs/math_function.hpp"

class Math_Function;
class Matrix_Solver {

	//from math_libs
	Math_Function math_f;

	//for LU decomposition
	struct LU_Matrices{
		std::vector<std::vector<double>> lhs_matrix_L;
		std::vector<std::vector<double>> lhs_matrix_U;
		std::vector<double> y;
	}LU_matrices;

	//internal functions
	void matrix_solver_cut_matrix(Variables &vars, int max_node);
//	void matrix_solver_LU_decomposition(std::vector<std::vector<double>> lhs_matrix, LU_Matrices &LU_matrices, int max_node);
//	void matrix_solver_backward_substitution_L(LU_Matrices &LU_matrices, std::vector<double> rhs_matrix, int max_node);
//	void matrix_solver_backward_substitution_U(LU_Matrices LU_matrices, std::vector<double> &lhs_result, int max_node);
	std::vector<double> matrix_solver_solve_matrix(Variables &vars);

	public:
		void matrix_solver_main_computation(Variables &vars, int max_node);

};

#endif
