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
	std::cout << " ... RHS MATRIX: DONE\n";
	lhs_matrix	=	matrix_init_lhs_matrix_calc(airfoil_pars, pars, vars);
	std::cout << " ... LHS MATRIX: DONE\n";

	misc.print_to_file(airfoil_pars.N1, "n1.dat");
	misc.print_to_file(airfoil_pars.N2, "n2.dat");
	misc.print_to_file(airfoil_pars.x, "coord_x.dat");
	misc.print_to_file(airfoil_pars.y, "coord_y.dat");
	misc.print_to_file(airfoil_pars.s, "s.dat");
	misc.print_to_file(airfoil_pars.nx, "nx.dat");
	misc.print_to_file(airfoil_pars.ny, "ny.dat");
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

	//compute dphi_inf/dn for N-
	double dphi_inf_minus_dn	=	matrix_init_dphi_dn_calc(nx[0], ny[0], pars);	

	//compute dphi_inf/dn for N+
	double dphi_inf_plus_dn		=	matrix_init_dphi_dn_calc(nx[max_node - 1], ny[max_node - 1], pars);	
	
	//calculation of RHS
	std::vector<double> G_N_plus	=	g_calc.g1_calc(airfoil_pars, x[0], y[0], max_node);
	std::vector<double> G_N_minus	=	g_calc.g2_calc(airfoil_pars, x[0], y[0], max_node);

	//calculate sigma 
	std::vector<double> sigma_H_ij_phi_j	=	matrix_init_H_rhs_calc(airfoil_pars, pars);

	for (auto i = 0; i < max_node; i++) {
		rhs_matrix[i]	=	G_N_minus[i]*dphi_inf_minus_dn + G_N_plus[i]*dphi_inf_plus_dn - sigma_H_ij_phi_j[i];
	}

	return rhs_matrix;
}

//now compute another lhs (version two
std::vector<std::vector<double>> Matrix_Init::matrix_init_lhs_matrix_calc(Airfoil_Parameters airfoil_pars, Parameters pars, Variables vars) {
	
	//make local pars
	double max_node		=	pars.max_node;
	std::vector<std::vector<double>> lhs_matrix(max_node, std::vector<double> (max_node));//row x column
	std::vector<std::vector<double>> test_matrix(max_node, std::vector<double> (max_node));//row x column

	//make local airfoil_pars 
	std::vector<double> nx	=	airfoil_pars.nx;
	std::vector<double> ny	=	airfoil_pars.ny;
	std::vector<double> x	=	airfoil_pars.x;
	std::vector<double> y	=	airfoil_pars.y;
	std::vector<double> s	=	airfoil_pars.s;

	//calculate sigma 
	std::vector<double> sigma_H_ij	=	matrix_init_H_lhs_calc(airfoil_pars, pars);
	
	//make local grids for integration
	//sigma -> dummy variable for integration
	
	double sigma_i	=	-1.;
	double sigma_f	=	1.;
	int n		=	99999;
	std::vector<double> sigma(n);
	std::vector<double> f_sigma(n);
	for (auto i = 0; i < n; i++) {
		double zeta	=	static_cast<double>(i)/static_cast<double>(n-1);
		sigma[i]	=	sigma_i*(1-zeta) + sigma_f*zeta;
	}

	for (auto i = 0; i < max_node; i++) {
		for (auto j = 0; j < max_node; j++) {

//			if (j == 1) {
			if (j == max_node - 1) {
				
				lhs_matrix[i][j]	=	-1*sigma_H_ij[i];
				
			} else {
/*				//calculate g_ij_1
				double temp_g_ij_1;

				if (i == 0 || j == 0) {
					temp_g_ij_1	=	0; //it will be cut anyway, so 0 is the best deal
				} else if (j == max_node - 1 && i == 1) {
					temp_g_ij_1	=	(s[j-1]/(8*M_PI))*(1 - 2*log(s[j-1]));
				} else if (i == j + 1) {
					temp_g_ij_1	=	(s[j-1]/(8*M_PI))*(1 - 2*log(s[j-1]));
				} else if (i == j) {
					temp_g_ij_1	=	(s[j-1]/(8*M_PI))*(3 - 2*log(s[j-1]));
				} else {
					temp_g_ij_1	=	g_calc.g_ij_1_calc(s[j-1], x[i], x[j], x[j-1], y[i], y[j], y[j-1]);
				}

				//calculate g_ij_2
				double temp_g_ij_2;

				if (i == 0 || j == 0) {
					temp_g_ij_2	=	0; //it will be cut anyway, so 0 is the best deal
				} else if (i == 1 && j == 1) {	//g_11_2
					temp_g_ij_2	=	(s[max_node-2]/(8*M_PI))*(3 - 2*log(s[max_node-2]));
				} else if (j == 1 && i == max_node - 1) {	//g_1N_2
					temp_g_ij_2	=	(s[max_node-2]/(8*M_PI))*(3 - 2*log(s[max_node-2]));
				} else if (i == j - 1) {
					temp_g_ij_2	=	(s[j-2]/(8*M_PI))*(3 - 2*log(s[j-2]));
				} else if (i == j) {
					temp_g_ij_2	=	(s[j-2]/(8*M_PI))*(3 - 2*log(s[j-2]));
				} else if (j == 1){
					temp_g_ij_2	=	g_calc.g_ij_2_calc(s[max_node-1], x[i], x[max_node-1], x[max_node-2], y[i], y[max_node-1], y[max_node-2]);
				} else {
					temp_g_ij_2	=	g_calc.g_ij_2_calc(s[j-2], x[i], x[j-1], x[j-2], y[i], y[j-1], y[j-2]);
				}
*/
				//PROCESS FOR INTEGRAL
				//g_1
				for (auto k = 0; k < n; k++) {
					f_sigma[k] = g_calc.g_1_integration_function(sigma[k], x[i], x[j], x[j+1], y[i], y[j], y[j+1]);
				}
				misc.neumann_bc(f_sigma);
				//integration process
				double g_ij_1_integrated	=	-1*(s[j-1]/(8*M_PI))*math_f.integral_simpson(sigma, f_sigma);
				
				//g_2
				for (auto k = 0; k < n; k++) {
					f_sigma[k] = g_calc.g_2_integration_function(sigma[k], x[i], x[j-1], x[j], y[i], y[j-1], y[j]);
				}
				misc.neumann_bc(f_sigma);
					//integration process
				double g_ij_2_integrated	=	-1*(s[j-2]/(8*M_PI))*math_f.integral_simpson(sigma, f_sigma);

				lhs_matrix[i][j]	=	g_ij_1_integrated + g_ij_2_integrated;
				//lhs_matrix[i][j]	=	temp_g_ij_2	+	temp_g_ij_1;
			}
		}
	}

	misc.print_to_file(test_matrix, "test_matrix_1.dat");
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

//function for calculate Hij of RHS
std::vector<double> Matrix_Init::matrix_init_H_rhs_calc(Airfoil_Parameters airfoil_pars, Parameters pars) {

	//make local vars
	int max_node			=	pars.max_node;
	std::vector<double> x		=	airfoil_pars.x;
	std::vector<double> y		=	airfoil_pars.y;
	std::vector<double> s		=	airfoil_pars.s;
	std::vector<double> beta	=	airfoil_pars.beta;
	std::vector<double> nx		=	airfoil_pars.nx;
	std::vector<double> ny		=	airfoil_pars.ny;

	double angle_of_attack		=	pars.angle_of_attack;
	double v_inf			=	pars.v_freestream;
	
	//make local grids for integration
	//sigma -> dummy variable for integration
	
	double sigma_i	=	-1.;
	double sigma_f	=	1.;
	int n		=	99999;
	std::vector<double> sigma(n);
	std::vector<double> f_sigma(n);
	for (auto i = 0; i < n; i++) {
		double zeta	=	static_cast<double>(i)/static_cast<double>(n-1);
		sigma[i]	=	sigma_i*(1-zeta) + sigma_f*zeta;
	}

	std::vector<double> H(max_node);
	for (auto i = 0; i < max_node; i++) {
		double sum_H_phi;
		for (auto j = 0; j < max_node; j++) {

			//================================
			//======= calculate h_ij_1 =======
			//================================
			double temp_h_ij_1;
			if (i == j || j == 1) {
				temp_h_ij_1	=	0;
			} else {
				for (auto k = 0; k < n; k++) {
					f_sigma[k] = h_calc.h_1_integration_function(sigma[k], nx[j], x[i], x[j], x[j+1], ny[j], y[i], y[j], y[j+1]);
				}
				misc.neumann_bc(f_sigma);
				//integration process
				temp_h_ij_1	=	-1*(s[j-1]/(8*M_PI))*math_f.integral_simpson(sigma, f_sigma);
			}
			//================================
			//======= calculate h_ij_2 =======
			//================================
			double temp_h_ij_2;
			if (i == j || j == 1) {
			temp_h_ij_2	=	0;
			} else if (j == 1) {
//			if (j == 1) {
				for (auto k = 0; k < n; k++) {
					f_sigma[k] = h_calc.h_2_integration_function(sigma[k], nx[max_node-1], x[i], x[max_node-1], x[1], ny[max_node-1], y[i], y[max_node-1], y[1]);
				}
				misc.neumann_bc(f_sigma);
				//integration process
				temp_h_ij_2	=	-1*(s[max_node-2]/(8*M_PI))*math_f.integral_simpson(sigma, f_sigma);
		
			} else {

				for (auto k = 0; k < n; k++) {
					f_sigma[k] = h_calc.h_2_integration_function(sigma[k], nx[j-1], x[i], x[j-1], x[j], ny[j-1], y[i], y[j-1], y[j]);
				}
				misc.neumann_bc(f_sigma);
				//integration process
				temp_h_ij_2	=	-1*(s[j-2]/(8*M_PI))*math_f.integral_simpson(sigma, f_sigma);
			}

			//=================================
			//======= calculate c_delta =======
			//=================================
			double c_delta;
			if (j == i) {
				c_delta	=	beta[j]/(2*M_PI);
			} else {
				c_delta	=	0;
			}

			//calculate phi_j
			double phi_j	=	matrix_init_phi_calc(x[j], y[j], pars);	
//			double temp_sum	=	(temp_h_ij_1 + temp_h_ij_2 + c_delta)*phi_j;
			double temp_sum	=	(temp_h_ij_1 + temp_h_ij_2 + c_delta)*(y[j]*cos(angle_of_attack) - x[j]*sin(angle_of_attack));

			//sum
			sum_H_phi	=	sum_H_phi + temp_sum;
		}
		H[i]	=	-1*v_inf*sum_H_phi;
	}
	return H;
}

//function for calculate Hij of LHS
std::vector<double> Matrix_Init::matrix_init_H_lhs_calc(Airfoil_Parameters airfoil_pars, Parameters pars) {

	//calculate integral method
	//make local vars
	int max_node			=	pars.max_node;
	std::vector<double> x		=	airfoil_pars.x;
	std::vector<double> y		=	airfoil_pars.y;
	std::vector<double> s		=	airfoil_pars.s;
	std::vector<double> beta	=	airfoil_pars.beta;
	std::vector<double> nx		=	airfoil_pars.nx;
	std::vector<double> ny		=	airfoil_pars.ny;
	
	//make local grids for integration
	//sigma -> dummy variable for integration
	
	double sigma_i	=	-1.;
	double sigma_f	=	1.;
	int n		=	99999;
	std::vector<double> sigma(n);
	std::vector<double> f_sigma(n);
	for (auto i = 0; i < n; i++) {
		double zeta	=	static_cast<double>(i)/static_cast<double>(n-1);
		sigma[i]	=	sigma_i*(1-zeta) + sigma_f*zeta;
	}

	std::vector<double> H(max_node);
	for (auto i = 0; i < max_node; i++) {
		double sum_H_phi;
		for (auto j = 0; j < max_node; j++) {

			//================================
			//======= calculate h_ij_1 =======
			//================================
			double temp_h_ij_1;
			double temp_integration_h1;
			if (i == j || j == 1) {
				temp_h_ij_1	=	0;
			} else {
				for (auto k = 0; k < n; k++) {
					f_sigma[k] = h_calc.h_1_integration_function(sigma[k], nx[j], x[i], x[j], x[j+1], ny[j], y[i], y[j], y[j+1]);
				}
				misc.neumann_bc(f_sigma);

				//integration process
				temp_integration_h1	=	math_f.integral_simpson(sigma, f_sigma);
				temp_h_ij_1		=	-1*(s[j-1]/(8*M_PI))*temp_integration_h1;
			}
/*			if (j < 5 && j > 0 && i < 5 && i > 0) {
				std::cout << i << " " << j << " " << temp_h_ij_1 << " " << temp_integration_h1 << " " << x[i] << " " << y[i] << std::endl;
				if (j == 1 && i == 0) {
					misc.print_to_file(f_sigma, "f_sigma_0.dat");
				} else if (i == 1) {
					misc.print_to_file(f_sigma, "f_sigma_1.dat");
				} else if (i == 2) {
					misc.print_to_file(f_sigma, "f_sigma_2.dat");
				} else if (i == 3) {
					misc.print_to_file(f_sigma, "f_sigma_3.dat");
				} else if (i == 4) {
					misc.print_to_file(f_sigma, "f_sigma_4.dat");
				}
			}
*/			
			//================================
			//======= calculate h_ij_2 =======
			//================================
			double temp_h_ij_2;
			double temp_integration_h2;
			if (i == j || j == 1) {
				temp_h_ij_2	=	0;
			} else if (j == 1) {
				for (auto k = 0; k < n; k++) {
					f_sigma[k] = h_calc.h_2_integration_function(sigma[k], nx[max_node-1], x[i], x[max_node-1], x[1], ny[max_node-1], y[i], y[max_node-1], y[1]);
				}
				misc.neumann_bc(f_sigma);

				//integration process
				temp_integration_h2	=	math_f.integral_simpson(sigma, f_sigma);	
				temp_h_ij_2		=	-1*(s[max_node-2]/(8*M_PI))*temp_integration_h2;
		
			} else {

				for (auto k = 0; k < n; k++) {
					f_sigma[k] = h_calc.h_2_integration_function(sigma[k], nx[j-1], x[i], x[j-1], x[j], ny[j-1], y[i], y[j-1], y[j]);
				}
				misc.neumann_bc(f_sigma);

				//integration process
				temp_integration_h2	=	math_f.integral_simpson(sigma, f_sigma);	
				temp_h_ij_2	=	-1*(s[j-2]/(8*M_PI))*temp_integration_h2;
			}

			//=================================
			//======= calculate c_delta =======
			//=================================
			double c_delta;
			if (j == i) {
				c_delta	=	beta[j]/(2*M_PI);
			} else {
				c_delta	=	0;
			}

			//calculate phi_j
			double phi_j	=	matrix_init_phi_calc(x[j], y[j], pars);	
			double temp_sum	=	(temp_h_ij_1 + temp_h_ij_2 + c_delta);
			

			//sum
			sum_H_phi	=	sum_H_phi + temp_sum;
		}
		H[i]	=	sum_H_phi;
	}
	return H;
}

//BELOW IS THE SECOND VERSION OF g1 and g2
//first: a_calc
double Matrix_Init::G_Calc::a_calc(double x_ref, double x_j, double x_j1, double y_ref, double y_j, double y_j1) {
	
	return 0.25*(pow(x_j,2) + pow(x_j1,2)) + 0.5*x_j1*x_ref + 0.25*(pow(y_j,2) + pow(y_j1,2)) + 0.5*y_j1*y_ref;
}

//second: b_calc
double Matrix_Init::G_Calc::b_calc(double x_ref, double x_j, double x_j1, double y_ref, double y_j, double y_j1) {
	
	double temp_1 = 0.5*(pow(x_j1,2) - pow(x_j,2)) - x_j1*x_ref + x_j*x_ref;
	double temp_2 = 0.5*(pow(y_j1,2) - pow(y_j,2)) - y_j1*y_ref + y_j*y_ref;

	return temp_1 + temp_2;
}

//second: c_calc
double Matrix_Init::G_Calc::c_calc(double x_ref, double x_j, double x_j1, double y_ref, double y_j, double y_j1) {
	
	double temp_1 = 0.25*(pow(x_j1,2) + pow(x_j,2)) + pow(x_ref,2) + 0.5*x_j1*x_ref - x_j1*x_ref - x_j*x_ref;
	double temp_2 = 0.25*(pow(y_j1,2) + pow(y_j,2)) + pow(y_ref,2) + 0.5*y_j1*y_ref - y_j1*y_ref - y_j*y_ref;

	return temp_1 + temp_2;
}

//compute I1
double Matrix_Init::G_Calc::I_1_calc(double a, double b, double c) {
	
	double delta 	= 	math_f.discriminant(a, b,c);
//	if (delta > 0) {
//		delta = delta*(-1);
//	
//	} else if (delta == 0) {
//		delta = delta - 0.0000001;
//	}

	double temp_1 	=	2/sqrt(-1*delta); 
	double temp_2	=	atan((sqrt(-1*delta))/(c - a));

	return temp_1*temp_2;
}

//compute I2
double Matrix_Init::G_Calc::I_2_calc(double a, double b, double c) {

	double delta 	= 	math_f.discriminant(a, b,c);
//	if (delta > 0) {
//		delta = delta*(-1);
//	} else if (delta == 0) {
//		delta = delta - 0.0000001;
//	}
	double temp_1	=	1/(2*a);
	double temp_2;
	if (a + b + c == 0) {
//		a	=	a + 0.00000001;
		temp_2	=	log((a + b + c)/(a - b + c));
	} else {
		temp_2	=	log((a + b + c)/(a - b + c));
	}
	double temp_3	=	b/(a*sqrt(-1*delta));
	double temp_4	=	atan((sqrt(-1*delta))/(c - a));
	

	return (temp_1*temp_2) - (temp_3*temp_4); 
}

//compute I3
double Matrix_Init::G_Calc::I_3_calc(double a, double b, double c, double I_1, double I_2) {

	if (a + b + c == 0) {
//		a	=	a + 0.00000001;
		//return temp - 4 + b*I_2 + 2*c*I_1;
		return log((a + b + c)*(a - b + c)) - 4 + b*I_2 + 2*c*I_1; 
	} else {
		return log((a + b + c)*(a - b + c)) - 4 + b*I_2 + 2*c*I_1; 
	}
}

//compute I4
double Matrix_Init::G_Calc::I_4_calc(double a, double b, double c, double I_3) {
	
	double temp_1	=	1/(2*a);
	double temp_2;
	if (a + b + c == 0) {
//		a	=	a + 0.00001;
		temp_2	=	log((a + b + c)/(a - b + c));
	} else {
		temp_2	=	log((a + b + c)/(a - b + c));
	}
	double temp_3	=	(a - b + c)*log(a - b + c);
	double temp_4	=	2*b;
	double temp_5	=	(b*I_3)/(2*a);

	return temp_1*(temp_2 - temp_3 - temp_4) - temp_5;
}

//compute general formula of g_ij_1
double Matrix_Init::G_Calc::g_ij_1_calc(double s, double x_ref, double x_j, double x_j1, double y_ref, double y_j, double y_j1) {

	//calculate a b and c
	double a	=	a_calc(x_ref, x_j, x_j1, y_ref, y_j, y_j1); 
	double b	=	b_calc(x_ref, x_j, x_j1, y_ref, y_j, y_j1); 
	double c	=	c_calc(x_ref, x_j, x_j1, y_ref, y_j, y_j1); 

	//compute I_1 I_2 I_3 and I_4
	double I_1	=	I_1_calc(a, b, c);
	double I_2	=	I_2_calc(a, b, c);
	double I_3	=	I_3_calc(a, b, c, I_1, I_2);
	double I_4	=	I_4_calc(a, b, c, I_3);

	//compute g_ij
	
	double temp_1	=	-1*s/(8*M_PI);
	double temp_2	=	0.5*(I_3 - I_4);

	return temp_1*temp_2;
}

//compute general formula of g_ij_2
double Matrix_Init::G_Calc::g_ij_2_calc(double s, double x_ref, double x_j, double x_j1, double y_ref, double y_j, double y_j1) {

	//calculate a b and c
	double a	=	a_calc(x_ref, x_j, x_j1, y_ref, y_j, y_j1); 
	double b	=	b_calc(x_ref, x_j, x_j1, y_ref, y_j, y_j1); 
	double c	=	c_calc(x_ref, x_j, x_j1, y_ref, y_j, y_j1); 

	//compute I_1 I_2 I_3 and I_4
	double I_1	=	I_1_calc(a, b, c);
	double I_2	=	I_2_calc(a, b, c);
	double I_3	=	I_3_calc(a, b, c, I_1, I_2);
	double I_4	=	I_4_calc(a, b, c, I_3);

	//compute g_ij
	
	double temp_1	=	-1*s/(8*M_PI);
	double temp_2	=	0.5*(I_3 + I_4);

	return temp_1*temp_2;
}
//g1_calc for G_Calc nested class
std::vector<double> Matrix_Init::G_Calc::g1_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node) {

	std::vector<double> g1(max_node);

	//make local variables
	std::vector<double> &s	=	airfoil_pars.s;
	std::vector<double> &x	=	airfoil_pars.x;
	std::vector<double> &y	=	airfoil_pars.y;

	//non-first-edge first
	for (auto i = 1; i < max_node; i++) {
		g1[i]	=	g_ij_1_calc(s[i-1], x_ref, x[i], x[i-1], y_ref, y[i], y[i-1]);
	}

	//0th node
	g1[0]	=	g1[max_node - 1];

	return g1;
}

//g2_calc for G_Calc nested class
std::vector<double> Matrix_Init::G_Calc::g2_calc(Airfoil_Parameters airfoil_pars, double x_ref, double y_ref, int max_node) {

	std::vector<double> g2(max_node);

	//make local variables
	std::vector<double> &s	=	airfoil_pars.s;
	std::vector<double> &x	=	airfoil_pars.x;
	std::vector<double> &y	=	airfoil_pars.y;

	//non-first-edge first
	for (auto i = 1; i < max_node; i++) {
		g2[i]	=	g_ij_2_calc(s[i-1], x_ref, x[i], x[i-1], y_ref, y[i], y[i-1]);
	}

	//0th node
	g2[0]	=	g2[max_node - 1];

	return g2;
}

//g_1_integration_function
double Matrix_Init::G_Calc::g_1_integration_function(double sigma, double x_ref, double x_j, double x_j1, double y_ref, double y_j, double y_j1) {
	
	double temp_1		=	1 - sigma;
	double temp_2		=	((0.5*(1 + sigma)*x_j1) + (0.5*(1 - sigma)*x_j) - (x_ref));
	double temp_3		=	((0.5*(1 + sigma)*y_j1) + (0.5*(1 - sigma)*y_j) - (y_ref));

	double r_ij		=	sqrt(pow(temp_2, 2) + pow(temp_3, 2));

	return temp_1*(log(r_ij));
}

//g_2_integration_function
double Matrix_Init::G_Calc::g_2_integration_function(double sigma, double x_ref, double x_j, double x_j1, double y_ref, double y_j, double y_j1) {
	
	double temp_1		=	1 + sigma;
	double temp_2		=	((0.5*(1 + sigma)*x_j1) + (0.5*(1 - sigma)*x_j) - (x_ref));
	double temp_3		=	((0.5*(1 + sigma)*y_j1) + (0.5*(1 - sigma)*y_j) - (y_ref));

	double r_ij		=	sqrt(pow(temp_2, 2) + pow(temp_3, 2));

	return temp_1*(log(r_ij));
}

//FOR H_CALC NESTED CLASS
//define integration function for h_1
double Matrix_Init::H_Calc::h_1_integration_function(double sigma, double nx, double x_ref, double x_j, double x_j1, double ny, double y_ref, double y_j, double y_j1) {
	double temp_1		=	1 - sigma;
	double temp_2		=	((0.5*(1 + sigma)*x_j1) + (0.5*(1 - sigma)*x_j) - (x_ref));
	double temp_3		=	((0.5*(1 + sigma)*y_j1) + (0.5*(1 - sigma)*y_j) - (y_ref));

	double numerator	=	temp_2*nx + temp_3*ny;
	double denominator	=	pow(temp_2, 2) + pow(temp_3, 2);

	return	temp_1*numerator/denominator;
}

//define integration function for h_2
double Matrix_Init::H_Calc::h_2_integration_function(double sigma, double nx, double x_ref, double x_j, double x_j1, double ny, double y_ref, double y_j, double y_j1) {
	double temp_1		=	1 + sigma;
	double temp_2		=	((0.5*(1 + sigma)*x_j1) + (0.5*(1 - sigma)*x_j) - (x_ref));
	double temp_3		=	((0.5*(1 + sigma)*y_j1) + (0.5*(1 - sigma)*y_j) - (y_ref));

	double numerator	=	temp_2*nx + temp_3*ny;
	double denominator	=	pow(temp_2, 2) + pow(temp_3, 2);

	return	temp_1*numerator/denominator;
}
