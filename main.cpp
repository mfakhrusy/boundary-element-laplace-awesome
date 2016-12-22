#include "global.hpp"
#include "pre_process/initialization.hpp"
#include "pre_process/airfoil.hpp"
#include "solver/matrix_init.hpp"
#include "solver/matrix_solver.hpp"
#include "misc/post_computation.hpp"

int main() {
	
	//define structs and classes
	Airfoil_Parameters airfoil_pars;
	Parameters pars;
	Variables vars;

	//class for solver
	Initialization inits;
	Airfoil airfoil;
	Matrix_Init matrix_init;
	Matrix_Solver matrix_solver;
	Post_Computation post_computation;

	//solving process
	inits.initialization_read_input(pars);
	airfoil.airfoil_main_computation(airfoil_pars, pars);
	matrix_init.matrix_init_main_computation(airfoil_pars, pars, vars);
	matrix_solver.matrix_solver_main_computation(vars);

	//post computation
	post_computation.post_main_calculation(airfoil_pars, pars, vars);
}
