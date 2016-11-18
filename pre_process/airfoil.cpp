#include "../global.hpp"
#include "airfoil.hpp"

void Airfoil::read_airfoil() {
	
	std::string airfoil_type;
	std::cout << "--------- Start Reading Airfoil -----------" << std::endl;
	std::cout << "Please input NACA AIRFOIL (e.g 0015 0012): ";
	std::cin >> airfoil_type;
	std::cout << std::endl << "Checking from databases..." << std::endl;

	// check airfoil exist or not in the database
	if (check_airfoil_from_databases(airfoil_type)) {
		std::cout << "OK, airfoil data has been read successfully..." << std::endl;
		std::string airfoil_title, airfoil_coor;
		//read airfoil data from the databases
		std::string airfoil_path = "airfoil_databases/" + airfoil_type;
		std::ifstream airfoil_input(airfoil_path.c_str()); //in C++11, no need .c_str()
		
		if (airfoil_input.is_open()) {
			//------- READ AIRFOIL DATA AND APPEND TO ARRAY -------
			std::string airfoil_title, airfoil_coordinate;
			std::getline(airfoil_input, airfoil_title); // read the first line (title!)
			int i = 0;
			while (std::getline(airfoil_input, airfoil_coordinate)) {
				std::istringstream iss(airfoil_coordinate);
				double a, b;
				iss >> a >> b;
				x.push_back(a);
				y.push_back(b);
				i++;
			}	
			nodes = i - 1; //define total nodes
			// sharpen the airfoil's trailing edge, for kutta condition purpose.
			y[nodes] = 0;
			y[0] = 0;

		} else {
			std::cout << "somehoew, the airfoil data is not exist!" << std::endl;
		}

	} else {
		std::cout << "Sorry, airfoil data is not found." << std::endl;
	}
}

bool Airfoil::check_airfoil_from_databases(std::string airfoil_input) {

	std::ifstream airfoil_databases;
	airfoil_databases.open("airfoil_databases/list_airfoil.dat");
	if (airfoil_databases.is_open()) {
		while (!airfoil_databases.eof()) {
			std::string read_per_line;
			airfoil_databases >> read_per_line;

			if (airfoil_input == read_per_line) {
				return true;
			} 
		}
		return false;
	}
	else {
		std::cout << "There isn't any file named 'list_airfoil.dat'" << std::endl;
		return false;
	}
}
