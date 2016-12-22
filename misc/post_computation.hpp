#ifndef POST_COMPUTATION_HPP
#define POST_COMPUTATION_HPP

class Miscellaneous;
class Post_Computation {

	std::vector<double> post_compute_v_tangential(Airfoil_Parameters airfoil_pars, Parameters pars, std::vector<double> lhs_result);
	std::vector<double> post_compute_cp(std::vector<double> v_tangential, double v_freestream);
	public:
		void post_main_calculation(Airfoil_Parameters airfoil_pars, Parameters pars, Variables &vars);

};

#endif
