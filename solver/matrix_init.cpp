#include "../global.hpp"
#include "matrix_init.hpp"

void Matrix_Init::matrix_init_main_computation(Airfoil_Parameters airfoil_pars, Parameters pars, Variables &vars) {

	//make local parameters
	double v_freestream			=	pars.v_freestream;
	double angle_of_attack			=	pars.angle_of_attack;
	int max_node				=	pars.max_node;

	std::vector<double> x			=	airfoil_pars.x;
	std::vector<double> y			=	airfoil_pars.y;

	//make local variables
	std::vector<double> &G_N_plus		=	vars.G_N_plus;
	std::vector<double> &G_N_minus		=	vars.G_N_minus;

	//calculation of RHS
	G_N_plus	=	matrix_init_g1_calc(airfoil_pars, x[0], y[0], max_node);
	G_N_minus	=	matrix_init_g2_calc(airfoil_pars, x[0], y[0], max_node);

}

//calculate (only one value) derivation of potential wrt normal (dphi/dn)
double Matrix_Init::matrix_init_calc(double nx, double ny, Parameters pars) {

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
	std::vector<double> &N2	=	airfoil_pars.N2;
	std::vector<double> &s	=	airfoil_pars.s;
	std::vector<double> &x	=	airfoil_pars.x;
	std::vector<double> &y	=	airfoil_pars.y;

	//non-first-edge first
	for (auto i = 1; i < max_node; i++) {

		g2[i]	=	-1*N2[i]*s[i-1]*math_f.green_f(x[i], x_ref, y[i], y_ref);
	}

	//0th node
	g2[0]	=	g2[max_node - 1];

	return g2;
}

