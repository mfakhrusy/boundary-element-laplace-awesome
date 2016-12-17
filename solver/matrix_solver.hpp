#ifndef MATRIX_SOLVER_HPP
#define MATRIX_SOLVER_HPP
#include "../math_libs/math_function.hpp"
#include "../misc/miscellaneous.hpp"

class Math_Function;
class Miscellaneous;
class Matrix_Solver {

	//friend classes
	Math_Function math_f;
	Miscellaneous misc;

	//internal functions
	void matrix_solver_cut_matrix(Variables &vars);
	std::vector<double> matrix_solver_solve_matrix(Variables &vars);

	public:
		void matrix_solver_main_computation(Variables &vars);

};

#endif
