#ifndef MATRIX_SOLVER_HPP
#define MATRIX_SOLVER_HPP
#include "../math_libs/math_function.hpp"

class Math_Function;
class Matrix_Solver {

	Math_Function math_f;
	void matrix_solver_cut_matrix(Variables &vars, int max_node);
	void matrix_solver_gauss_elimination(Variables &vars, int max_node);

	public:
		void matrix_solver_main_computation(Variables &vars, int max_node);

};

#endif
