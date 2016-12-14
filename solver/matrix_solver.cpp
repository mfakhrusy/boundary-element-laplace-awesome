#include "../global.hpp"
#include "matrix_solver.hpp"
	
void Matrix_Solver::matrix_solver_main_computation(Variables &vars, int max_node) {

	//cut the matrix
	matrix_solver_cut_matrix(vars, max_node);

}

//first, cut the matrix from 131x131 to 130x130 (lhs) and 131 to 130 (rhs)
void Matrix_Solver::matrix_solver_cut_matrix(Variables &vars, int max_node) {

	//make local vars
	std::vector<double> &rhs_matrix			=	vars.rhs_matrix;
	std::vector<std::vector<double>> &lhs_matrix	=	vars.lhs_matrix;

	//delete the rhs_first element
	std::ofstream output_file;
	output_file.open("output-temp/rhs_before.dat");
	for (auto i = 0; i < rhs_matrix.size(); i++) {
		output_file << i << " " << rhs_matrix[i] << std::endl;
	}
	output_file.close();
	
	//cut function start
	rhs_matrix.erase(rhs_matrix.begin(), rhs_matrix.begin() + 1); //delete the 0th node
	//cut function end 

	output_file.open("output-temp/rhs_after.dat");
	for (auto i = 0; i < rhs_matrix.size(); i++) {
		output_file << i << " " << rhs_matrix[i] << std::endl;
	}
	output_file.close();

	//delete the lhs first row and column

	output_file.open("output-temp/lhs_before.dat");
	for (auto i = 0; i < lhs_matrix.size(); i++) {
		for (auto j = 0; j < lhs_matrix[i].size(); j++) {
			output_file << lhs_matrix[i][j] << " ";
		}
		output_file << std::endl;
	}
	output_file.close();

	//cut function start
	lhs_matrix.erase(lhs_matrix.begin(), lhs_matrix.begin() + 1);
	for (auto i = 0; i < lhs_matrix.size(); i++) {
		lhs_matrix[i].erase(lhs_matrix[i].begin(), lhs_matrix[i].begin() + 1);
	}
	//cut function end 

	output_file.open("output-temp/lhs_after.dat");
	for (auto i = 0; i < lhs_matrix.size(); i++) {
		for (auto j = 0; j < lhs_matrix[i].size(); j++) {
			output_file << lhs_matrix[i][j] << " ";
		}
		output_file << std::endl;
	}
	output_file.close();

}
