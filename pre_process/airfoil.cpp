#include "../global.hpp"
#include "airfoil.hpp"

void Airfoil::read_airfoil() {
	
	std::string airfoil_type;
	std::cout << "--------- Start Reading Airfoil -----------" << std::endl;
	std::cout << "Please input NACA AIRFOIL (e.g 0015 0012): ";
	std::cin >> airfoil_type;
	std::cout << std::endl << "Checking from databases..." << std::endl;
	check_airfoil_from_databases(airfoil_type);
	std::cout << "OK, airfoil data has been read successfully..." << std::endl;
}

void Airfoil::check_airfoil_from_databases(std::string airfoil_input) {

	std::ifstream airfoil_databases;
	airfoil_databases.open("airfoil_databases/list_airfoil.dat");
	if (airfoil_databases.is_open()) {
		while (!airfoil_databases.eof()) {
			std::string read_per_line;
			airfoil_databases >> read_per_line;

			if (airfoil_input == read_per_line) {
				std::cout << "OK! Airfoil Exists, continue" << std::endl;
				break;
			} 

			if (airfoil_databases.eof()) {
				std::cout << "Airfoil is not exist in the database, sorry!" << std::endl;
			}
		}
		
	}
	else {
		std::cout << "There isn't any file named 'list_airfoil.dat'" << std::endl;
	}

}
