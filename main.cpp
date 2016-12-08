#include "global.hpp"
#include "pre_process/airfoil.hpp"

int main() {
	
	//define structs and classes
	Airfoil_Parameters airfoil_pars;
	Parameters pars;
	Airfoil airfoil;

	airfoil.airfoil_main_computation(airfoil_pars, pars);
}
