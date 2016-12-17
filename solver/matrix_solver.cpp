#include "../global.hpp"
#include "matrix_solver.hpp"
	
void Matrix_Solver::matrix_solver_main_computation(Variables &vars) {

	//make local vars
	std::vector<double> &lhs_result	=	vars.lhs_result;

	//cut the matrix
	matrix_solver_cut_matrix(vars);

	//solve it
	lhs_result	=	matrix_solver_solve_matrix(vars);
}

//first, cut the matrix from 131x131 to 130x130 (lhs) and 131 to 130 (rhs)
void Matrix_Solver::matrix_solver_cut_matrix(Variables &vars) {

	//make local vars
	std::vector<double> &rhs_matrix			=	vars.rhs_matrix;
	std::vector<std::vector<double>> &lhs_matrix	=	vars.lhs_matrix;

	const int max_node	=	rhs_matrix.size();

	//delete the rhs_first element
	misc.print_to_file(rhs_matrix, "rhs_before_cut.dat");
	
	//cut function start
	rhs_matrix.erase(rhs_matrix.begin(), rhs_matrix.begin() + 1); //delete the 0th node
	//cut function end 

	misc.print_to_file(rhs_matrix, "rhs_after_cut.dat");

	//delete the lhs first row and column

	misc.print_to_file(lhs_matrix, "lhs_before_cut.dat");

	//cut function start
	lhs_matrix.erase(lhs_matrix.begin(), lhs_matrix.begin() + 1);
	for (auto i = 0; i < lhs_matrix.size(); i++) {
		lhs_matrix[i].erase(lhs_matrix[i].begin(), lhs_matrix[i].begin() + 1);
	}
	//cut function end 

	misc.print_to_file(lhs_matrix, "lhs_after_cut.dat");

}

//gauss elimination part, make upper triangle.
std::vector<double> Matrix_Solver::matrix_solver_solve_matrix(Variables &vars) {
	//from math_libs
	
	//make local vars
	std::vector<std::vector<double>> lhs	=	vars.lhs_matrix;
	std::vector<double> rhs			=	vars.rhs_matrix;

	const int max_node	=	rhs.size();

	std::vector<double> lhs_result(max_node);

	//first combine lhs and rhs into 1 matrix;
	std::vector<std::vector<double>> lhs_temp	=	lhs;	//first, make a copy of lhs
	
	//iterating the value for combining the rhs_matrix into lhs_matrix
	for (auto i = 0; i < max_node; i++) {
		lhs_temp[i].push_back(rhs[i]);
	}

	//print the newly born matrix
	misc.print_to_file(lhs_temp, "lhs_plus_rhs.dat");

	for (auto i = 0; i < max_node; i++) {
	
		//search for the maximum in this column
		double max_E1	=	fabs(lhs_temp[i][i]);
		int max_row	=	i;
		for (auto k = i + 1; k < max_node; k++) {
			if (fabs(lhs_temp[k][i]) > max_E1) {
			
			max_E1	=	fabs(lhs_temp[k][i]);
			max_row	=	k;
			}
		}

		//swap maximum row with current row (column by column)
		for (auto k = 1; k < max_node; k++) {
			double temp		=	lhs_temp[max_row][k];
			lhs_temp[max_row][k]	=	lhs_temp[i][k];
			lhs_temp[i][k]		=	temp;
		}

		//make all rows below this one 0 in current column
		for (auto k = i + 1; k < max_node; k++) {
			double temp	=	-1*lhs_temp[k][i]/lhs_temp[i][i];
			for (auto j = i; j < max_node + 1; j++) {
				if (i == j) {
					lhs_temp[k][j]	=	0;
				} else {
					lhs_temp[k][j]	+=	temp*lhs_temp[i][j];
				}
			}
		}

	}

	misc.print_to_file(lhs_temp, "lhs_after_gauss.dat");

	//solve equation Ax=b for an upper triangular matrix lhs_temp
	
	for (auto i = max_node - 1; i >= 0; i--) {
		lhs_result[i]	=	lhs_temp[i][max_node]/lhs_temp[i][i];
		for (auto j = i - 1; j >= 0; j--) {
			lhs_temp[j][max_node]	-=	lhs_temp[j][i]*lhs_result[i];
		}
	}

	misc.print_to_file(lhs_result, "lhs_result.dat");

	return lhs_result;
}
