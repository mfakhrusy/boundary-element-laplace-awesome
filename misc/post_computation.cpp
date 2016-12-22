#include "../global.hpp"
#include "post_computation.hpp"
#include "../misc/miscellaneous.hpp"

void Post_Computation::post_main_calculation(Airfoil_Parameters airfoil_pars, Parameters pars, Variables &vars) {
	
	Miscellaneous misc;
	//make local vars
	std::vector<double> lhs_result		=	vars.lhs_result;
	std::vector<double> &v_tangential	=	vars.v_tangential;
	std::vector<double> &cp			=	vars.cp;
	const double v_freestream		=	pars.v_freestream;

	//compute v_tangential
	v_tangential	=	post_compute_v_tangential(airfoil_pars, pars, lhs_result);

	//compute cp distribution
	cp		=	post_compute_cp(v_tangential, v_freestream);


	//print to file
	misc.print_to_file(cp, "cp_distribution.dat");
	misc.print_to_file(v_tangential, "v_tangential_distribution.dat");
}

std::vector<double> Post_Computation::post_compute_v_tangential(Airfoil_Parameters airfoil_pars, Parameters pars, std::vector<double> lhs_result) {
	
	const int n	=	lhs_result.size();

	//make local parameters
	const double aoa	=	pars.angle_of_attack;
	const double v_inf	=	pars.v_freestream;

	std::vector<double> nx	=	airfoil_pars.nx;
	std::vector<double> ny	=	airfoil_pars.ny;
	
	//cut the 0th matrix -> because lhs_result contains 1 - 130 node, nx and ny contains 0 - 130 node.
	nx.erase(nx.begin(), nx.begin() + 1); //delete the 0th node
	ny.erase(ny.begin(), ny.begin() + 1); //delete the 0th node

	std::vector<double> v_tangential(n);
	for (auto i = 0; i < n; i++) {
		v_tangential[i]	=	lhs_result[i] + v_inf*(ny[i]*cos(aoa) - nx[i]*sin(aoa));
	}

	return v_tangential;
}

std::vector<double> Post_Computation::post_compute_cp(std::vector<double> v_tangential, const double v_freestream) {
	
	const int n	=	v_tangential.size();

	std::vector<double> cp(n);
	for (auto i = 0; i < n; i++) {
		cp[i]	=	1 - (pow(v_tangential[i], 2)/pow(v_freestream, 2));
	}
	return cp;
}
