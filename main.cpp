//the main reference to be used is paper by Iannelli et. al 1997: "A Kutta
// Condition Enforcing BEM Technology for Airfoil Aerodynamics

#include "global.hpp"
#include "pre_process/airfoil.hpp"

int main() {
	
	//inputs
	double v_stream = 10;
	double angle_of_attack = 0;

	//pre processing section + define local variable nodes, x, and y.
	Airfoil airfoil;
	airfoil.read_airfoil();
	double nodes = airfoil.nodes; // 0 index from upper TE, counter clockwise to lower TE.
	std::vector<double> x = airfoil.x;
	std::vector<double> y = airfoil.y;
	
	//solving process

}
