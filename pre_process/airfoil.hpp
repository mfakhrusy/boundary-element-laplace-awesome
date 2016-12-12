#ifndef AIRFOIL_HPP
#define AIRFOIL_HPP
#include "../math_libs/math_function.hpp"

class Math_Function;
class Airfoil {
	
	Math_Function math_f;

	//private parts -- functions
	void airfoil_read(std::vector<double> &x, std::vector<double> &y, int &max_node);	//max_node is initialized here
	bool airfoil_check_from_databases(std::string airfoil_input);
	std::vector<double> airfoil_element_length_calc(Airfoil_Parameters airfoil_pars, int max_node);
	std::vector<double> airfoil_beta_calc(Airfoil_Parameters airfoil_pars, int max_node);
	std::vector<double> airfoil_nx_calc(Airfoil_Parameters airfoil_pars, int max_node);
	std::vector<double> airfoil_ny_calc(Airfoil_Parameters airfoil_pars, int max_node);
	std::vector<double> airfoil_eta_calc(Airfoil_Parameters airfoil_pars, int max_node);
	std::vector<double> airfoil_N1_calc(Airfoil_Parameters airfoil_pars, int max_node);
	std::vector<double> airfoil_N2_calc(Airfoil_Parameters airfoil_pars, int max_node);
	

	// x and y for airfoil coordinates, read it from upper trailing edge, rotate clockwise until lower trailing edge

	public:
		void airfoil_main_computation(Airfoil_Parameters &airfoil_pars, Parameters &pars);

};
#endif
