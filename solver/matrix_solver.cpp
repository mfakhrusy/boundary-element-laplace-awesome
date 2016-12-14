#include "../global.hpp"
#include "matrix_solver.hpp"
	
void Matrix_Solver::matrix_solver_main_computation(Variables &vars, int max_node) {

	//make local vars
	std::vector<double> &lhs_result	=	vars.lhs_result;
	//cut the matrix
	matrix_solver_cut_matrix(vars, max_node);

	//solve it
	lhs_result	=	matrix_solver_solve_matrix(vars);
	
}

//first, cut the matrix from 131x131 to 130x130 (lhs) and 131 to 130 (rhs)
void Matrix_Solver::matrix_solver_cut_matrix(Variables &vars, int max_node) {

	//make local vars
	std::vector<double> &rhs_matrix			=	vars.rhs_matrix;
	std::vector<std::vector<double>> &lhs_matrix	=	vars.lhs_matrix;

	//delete the rhs_first element
	std::ofstream output_file;
	output_file.open("output-temp/rhs_before_cut.dat");
	for (auto i = 0; i < rhs_matrix.size(); i++) {
		output_file << i << " " << rhs_matrix[i] << std::endl;
	}
	output_file.close();
	
	//cut function start
	rhs_matrix.erase(rhs_matrix.begin(), rhs_matrix.begin() + 1); //delete the 0th node
	//cut function end 

	output_file.open("output-temp/rhs_after_cut.dat");
	for (auto i = 0; i < rhs_matrix.size(); i++) {
		output_file << i << " " << rhs_matrix[i] << std::endl;
	}
	output_file.close();

	//delete the lhs first row and column

	output_file.open("output-temp/lhs_before_cut.dat");
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

	output_file.open("output-temp/lhs_after_cut.dat");
	for (auto i = 0; i < lhs_matrix.size(); i++) {
		for (auto j = 0; j < lhs_matrix[i].size(); j++) {
			output_file << lhs_matrix[i][j] << " ";
		}
		output_file << std::endl;
	}
	output_file.close();

}

//gauss elimination part, make upper triangle.
std::vector<double> Matrix_Solver::matrix_solver_solve_matrix(Variables &vars) {
	
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
	std::ofstream output_file;
	output_file.open("output-temp/lhs_plus_rhs.dat");
	for (auto i = 0; i < lhs_temp.size(); i++) {
		for (auto j = 0; j < lhs_temp[i].size(); j++) {
			output_file << lhs_temp[i][j] << " ";
		}
		output_file << std::endl;
	}
	output_file.close();

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

	output_file.open("output-temp/lhs_after_gauss.dat");
	for (auto i = 0; i < lhs_temp.size(); i++) {
		for (auto j = 0; j < lhs_temp[i].size(); j++) {
			output_file << lhs_temp[i][j] << " ";
		}
		output_file << std::endl;
	}
	output_file.close();

	//solve equation Ax=b for an upper triangular matrix lhs_temp
	
	for (auto i = max_node - 1; i >= 0; i--) {
		lhs_result[i]	=	lhs_temp[i][max_node]/lhs_temp[i][i];
		for (auto j = i - 1; j >= 0; j--) {
			lhs_temp[j][max_node]	-=	lhs_temp[j][i]*lhs_result[i];
		}
	}

	output_file.open("output-temp/lhs_result.dat");
	for (auto i = 0; i < lhs_result.size(); i++) {
		output_file << i << " " << lhs_result[i];
		output_file << std::endl;
	}
	output_file.close();


	return lhs_result;
}
