#include "../global.hpp"
#include "initialization.hpp"

void Initialization::initialization_read_input(Parameters &pars) {

	//make local variables
	double &v_freestream	=	pars.v_freestream;
	double &angle_of_attack	=	pars.angle_of_attack;
	
	//read file input
	std::ifstream input_parameters_file;
	input_parameters_file.open("input/input.dat");

	double value;
	std::string value_title;
	int count = 0;

	while (input_parameters_file >> value >> value_title) {
		count++;

		switch(count) {
			case 1:
				v_freestream	=	value;
				break;
			case 2:
				angle_of_attack	=	value;
				break;
			default:
				std::cout << "Nothing to read!\n\n";
		}
	}
}
