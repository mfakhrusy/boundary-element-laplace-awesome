#include "../global.hpp"
#include "matrix_init.hpp"

void Matrix_Init::matrix_init_main_computation(Airfoil_Parameters airfoil_pars, Parameters pars, Variables &vars) {

	//make local parameters
	double v_freestream			=	pars.v_freestream;
	double angle_of_attack			=	pars.angle_of_attack;
	int max_node				=	pars.max_node;

	std::vector<double> x			=	airfoil_pars.x;
	std::vector<double> y			=	airfoil_pars.y;

	//make local vars
	std::vector<double> &rhs_matrix			=	vars.rhs_matrix;
	std::vector<std::vector<double>> &lhs_matrix	=	vars.lhs_matrix;

	//process
	rhs_matrix	=	matrix_init_rhs_matrix_calc(airfoil_pars, pars, vars);
	lhs_matrix	=	matrix_init_lhs_matrix_calc(airfoil_pars, pars, vars);
}

//main computation for RHS of equation
std::vector<double> Matrix_Init::matrix_init_rhs_matrix_calc(Airfoil_Parameters airfoil_pars, Parameters pars, Variables vars) {

	//make local pars
	double max_node		=	pars.max_node;
	std::vector<double> rhs_matrix(max_node);

	//make local airfoil_pars 
	std::vector<double> nx	=	airfoil_pars.nx;
	std::vector<double> ny	=	airfoil_pars.ny;
	std::vector<double> x	=	airfoil_pars.x;
	std::vector<double> y	=	airfoil_pars.y;

	//make local variables

	//compute dphi_inf/dn for N-
	double dphi_inf_minus_dn	=	matrix_init_dphi_dn_calc(nx[0], ny[0], pars);	

	//compute dphi_inf/dn for N+
	double dphi_inf_plus_dn		=	matrix_init_dphi_dn_calc(nx[max_node - 1], ny[max_node - 1], pars);	
	
	//calculation of RHS
	std::vector<double> G_N_plus	=	matrix_init_g1_calc(airfoil_pars, x[0], y[0], max_node);
	std::vector<double> G_N_minus	=	matrix_init_g2_calc(airfoil_pars, x[0], y[0], max_node);

	//calculate sigma 
	std::vector<double> sigma_H_ij_phi_j	=	matrix_init_H_rhs_calc(airfoil_pars, pars);

	for (auto i = 0; i < max_node; i++) {
		rhs_matrix[i]	=	G_N_minus[i]*dphi_inf_minus_dn + G_N_plus[i]*dphi_inf_plus_dn - sigma_H_ij_phi_j[i];
	}

	return rhs_matrix;
}

//main computation for LHS of equation
std::vector<std::vector<double>> Matrix_Init::matrix_init_lhs_matrix_calc(Airfoil_Parameters airfoil_pars, Parameters pars, Variables vars) {

	//make local pars
	double max_node		=	pars.max_node;
	std::vector<std::vector<double>> lhs_matrix(max_node, std::vector<double> (max_node));//row x column

	//make local airfoil_pars 
	std::vector<double> nx	=	airfoil_pars.nx;
	std::vector<double> ny	=	airfoil_pars.ny;
	std::vector<double> x	=	airfoil_pars.x;
	std::vector<double> y	=	airfoil_pars.y;

	//calculate sigma 
	std::vector<double> sigma_H_ij	=	matrix_init_H_lhs_calc(airfoil_pars, pars);

	for (auto i = 0; i < max_node; i++) {
		for (auto j = 0; j < max_node; j++) {

			if (j == max_node - 1) {
				
				lhs_matrix[i][j]	=	-1*sigma_H_ij[i];
				
			} else {
			
				//calculate g_ij_2
				std::vector<double> g_ij_2	=	matrix_init_g2_calc(airfoil_pars, x[j], y[j], max_node);
				double temp_g_ij_2		=	g_ij_2[i];

				//calculate g_i(j+1)_1
				std::vector<double> g_ij_1;
				double temp_g_ij_1;
				if (j != max_node - 1) {	//basically said: if j = max_node - 1 -> j = 1;
					g_ij_1	=	matrix_init_g1_calc(airfoil_pars, x[j+1], y[j+1], max_node);
					temp_g_ij_1			=	g_ij_1[i];
				} else {
					//impossible to reach LOLOLOLOLOLWKWKWK
					g_ij_1	=	matrix_init_g1_calc(airfoil_pars, x[1], y[1], max_node);
					temp_g_ij_1			=	g_ij_1[i];
				}

				if (j == 5 && i == 99) {
					for (auto k = 0; k < max_node; k++) {
						std::cout << k << " " << g_ij_1[k] << " " << g_ij_2[k] << std::endl;
					}
				}

				lhs_matrix[i][j]	=	temp_g_ij_2	+	temp_g_ij_1;

			}
		}
	}

	return lhs_matrix;
}

//calculate phi
double Matrix_Init::matrix_init_phi_calc(double x, double y, Parameters pars) {
	
	//make local vars
	double aoa	=	pars.angle_of_attack;
	double v_inf	=	pars.v_freestream;

	return v_inf*(y*cos(aoa) - x*sin(aoa));
}

//calculate (only one value) derivation of potential wrt normal (dphi/dn)
double Matrix_Init::matrix_init_dphi_dn_calc(double nx, double ny, Parameters pars) {

	double dphi_dn;

	//make local vars
	double aoa	=	pars.angle_of_attack;
	double v_inf	=	pars.v_freestream;

	dphi_dn		=	-1*v_inf*(nx*cos(aoa) - ny*sin(aoa));

	return dphi_dn;
}

//function for calculate g1 (in i value, need to loop for j value)
std::vector<double> Matrix_Init::matrix_init_g1_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node) {

	std::vector<double> g1(max_node);

	//make local variables
	std::vector<double> &N1	=	airfoil_pars.N1;
	std::vector<double> &s	=	airfoil_pars.s;
	std::vector<double> &x	=	airfoil_pars.x;
	std::vector<double> &y	=	airfoil_pars.y;

	//non-first-edge first
	for (auto i = 1; i < max_node; i++) {

		g1[i]	=	-1*N1[i]*s[i-1]*math_f.green_f(x[i], x_ref, y[i], y_ref);

	}

	//0th node
	g1[0]	=	g1[max_node - 1];

	return g1;
}

//function for calculate g2 (in i value, need to loop for j value)
std::vector<double> Matrix_Init::matrix_init_g2_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node) {

	std::vector<double> g2(max_node);

	//make local variables
	std::vector<double> N2	=	airfoil_pars.N2;
	std::vector<double> s	=	airfoil_pars.s;
	std::vector<double> x	=	airfoil_pars.x;
	std::vector<double> y	=	airfoil_pars.y;

	//non-first-edge first
	for (auto i = 1; i < max_node; i++) {

		g2[i]	=	-1*N2[i]*s[i-1]*math_f.green_f(x[i], x_ref, y[i], y_ref);
	}

	//0th node
	g2[0]	=	g2[max_node - 1];

	return g2;
}

//function for calculate h1(in i value, need to loop through for j value)
std::vector<double> Matrix_Init::matrix_init_h1_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node) {
	
	std::vector<double> h1(max_node);

	//make local vars
	std::vector<double> N1	=	airfoil_pars.N1;
	std::vector<double> s	=	airfoil_pars.s;
	std::vector<double> x	=	airfoil_pars.x;
	std::vector<double> y	=	airfoil_pars.y;
	std::vector<double> ny	=	airfoil_pars.ny;
	std::vector<double> nx	=	airfoil_pars.nx;

	//non-first-edge first
	for (auto i = 1; i < max_node; i++) {

		h1[i]	=	-1*N1[i]*s[i-1]*math_f.dgreen_f_dn(x[i], x_ref, nx[i], y[i], y_ref, ny[i]);
	}

	//0th node
	h1[0]	=	h1[max_node - 1];

	return h1;
}

//function for calculate h2(in i value, need to loop through for j value)
std::vector<double> Matrix_Init::matrix_init_h2_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node) {
	
	std::vector<double> h2(max_node);

	//make local vars
	std::vector<double> N2	=	airfoil_pars.N2;
	std::vector<double> s	=	airfoil_pars.s;
	std::vector<double> x	=	airfoil_pars.x;
	std::vector<double> y	=	airfoil_pars.y;
	std::vector<double> ny	=	airfoil_pars.ny;
	std::vector<double> nx	=	airfoil_pars.nx;

	//non-first-edge first
	for (auto i = 1; i < max_node; i++) {

		h2[i]	=	-1*N2[i]*s[i-1]*math_f.dgreen_f_dn(x[i], x_ref, nx[i], y[i], y_ref, ny[i]);
	}

	//0th node
	h2[0]	=	h2[max_node - 1];

	return h2;
}

//function for calculate Hij of RHS
std::vector<double> Matrix_Init::matrix_init_H_rhs_calc(Airfoil_Parameters airfoil_pars, Parameters pars) {

	//make local vars
	int max_node			=	pars.max_node;
	double v_freestream		=	pars.v_freestream;
	double angle_of_attack		=	pars.angle_of_attack;
	std::vector<double> x		=	airfoil_pars.x;
	std::vector<double> y		=	airfoil_pars.y;
	std::vector<double> s		=	airfoil_pars.s;
	std::vector<double> beta	=	airfoil_pars.beta;
	
	std::vector<double> H(max_node);
	for (auto i = 0; i < max_node; i++) {

		double sum_H_phi;
		for (auto j = 0; j < max_node; j++) {
			//calculate c_delta_ij
			double c_delta_ij;	
			if (j == i) {
				c_delta_ij	=	beta[i]/(2*M_PI);
			}
			else {
				c_delta_ij	=	0;
			}

			//calculate h_ij_2
			std::vector<double> h_ij_2	=	matrix_init_h2_calc(airfoil_pars, x[j], y[j], max_node);
			double temp_h_ij_2		=	h_ij_2[i];

			//calculate h_i(j+1)_1
			double temp_h_ij_1;
			if (j != max_node - 1) {
				std::vector<double> h_ij_1	=	matrix_init_h1_calc(airfoil_pars, x[j+1], y[j+1], max_node);
				temp_h_ij_1			=	h_ij_1[i];
			} else {
				std::vector<double> h_ij_1	=	matrix_init_h1_calc(airfoil_pars, x[1], y[1], max_node);
				temp_h_ij_1			=	h_ij_1[i];
			}

			//calculate phi_j
			double phi_j	=	matrix_init_phi_calc(x[j], y[j], pars);	
			
			//calculate the partial summation
			double temp_sum	=	(c_delta_ij - temp_h_ij_2 - temp_h_ij_1)*phi_j;

			//sum
			sum_H_phi	=	sum_H_phi + temp_sum;
		}
		
		H[i]	=	sum_H_phi;

	}
	return H;
}

//function for calculate Hij of LHS
std::vector<double> Matrix_Init::matrix_init_H_lhs_calc(Airfoil_Parameters airfoil_pars, Parameters pars) {

	//make local vars
	int max_node			=	pars.max_node;
	double v_freestream		=	pars.v_freestream;
	double angle_of_attack		=	pars.angle_of_attack;
	std::vector<double> x		=	airfoil_pars.x;
	std::vector<double> y		=	airfoil_pars.y;
	std::vector<double> s		=	airfoil_pars.s;
	std::vector<double> beta	=	airfoil_pars.beta;
	
	std::vector<double> H(max_node);
	for (auto i = 0; i < max_node; i++) {

		double sum_H_phi;
		for (auto j = 0; j < max_node; j++) {
			//calculate c_delta_ij
			double c_delta_ij;	
			if (j == i) {
				c_delta_ij	=	beta[i]/(2*M_PI);
			}
			else {
				c_delta_ij	=	0;
			}

			//calculate h_ij_2
			std::vector<double> h_ij_2	=	matrix_init_h2_calc(airfoil_pars, x[j], y[j], max_node);
			double temp_h_ij_2		=	h_ij_2[i];

			//calculate h_i(j+1)_1
			double temp_h_ij_1;
			if (j != max_node - 1) {
				std::vector<double> h_ij_1	=	matrix_init_h1_calc(airfoil_pars, x[j+1], y[j+1], max_node);
				temp_h_ij_1			=	h_ij_1[i];
			} else {
				std::vector<double> h_ij_1	=	matrix_init_h1_calc(airfoil_pars, x[1], y[1], max_node);
				temp_h_ij_1			=	h_ij_1[i];
			}

			//calculate the partial summation
			double temp_sum	=	(c_delta_ij - temp_h_ij_2 - temp_h_ij_1); //the difference from RHS is here, no phi_j

			//sum
			sum_H_phi	=	sum_H_phi + temp_sum;
		}
		
		H[i]	=	sum_H_phi;

	}
	return H;
}
