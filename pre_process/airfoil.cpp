#include "../global.hpp"
#include "airfoil.hpp"

void Airfoil::airfoil_main_computation(Airfoil_Parameters &airfoil_pars, Parameters &pars) {
	
	//make local variables
	int &max_node			=	pars.max_node;		//HAVE TO BE REFERENCE TO MAX_NODE BECAUSE AIRFOIL_READ IS INPUTING VALUE TO IT TNX
	std::vector<double> &x		=	airfoil_pars.x;
	std::vector<double> &y		=	airfoil_pars.y;
	std::vector<double> &s		=	airfoil_pars.s;
	std::vector<double> &beta	=	airfoil_pars.beta;
	std::vector<double> &nx		=	airfoil_pars.nx;
	std::vector<double> &ny		=	airfoil_pars.ny;
	std::vector<double> &eta	=	airfoil_pars.eta;
	std::vector<double> &N1		=	airfoil_pars.N1;
	std::vector<double> &N2		=	airfoil_pars.N2;

	//read the coordinate from databases
	airfoil_read(x, y, max_node);
	
	//calculate length
	s	=	airfoil_element_length_calc(airfoil_pars, max_node);

	//calculate beta
	beta	=	airfoil_beta_calc(airfoil_pars, max_node);

	//calculate nx and ny
	nx	=	airfoil_nx_calc(airfoil_pars, max_node);
	ny	=	airfoil_ny_calc(airfoil_pars, max_node);

	//calculate eta, N1, and N2
	eta	=	airfoil_eta_calc(airfoil_pars, max_node);
	N1	=	airfoil_N1_calc(airfoil_pars, max_node);
	N2	=	airfoil_N2_calc(airfoil_pars, max_node);
}

void Airfoil::airfoil_read(std::vector<double> &x, std::vector<double> &y, int &max_node) {

	std::string airfoil_type;
	std::cout << "--------- Start Reading Airfoil -----------" << std::endl;
	std::cout << "Please input NACA AIRFOIL (e.g 0015 0012): ";
	std::cin >> airfoil_type;
	std::cout << std::endl << "Checking from databases..." << std::endl;

	// check airfoil exist or not in the database
	if (airfoil_check_from_databases(airfoil_type)) {
		std::cout << "OK, airfoil data has been read successfully..." << std::endl;
		std::string airfoil_title, airfoil_coor;
		//read airfoil data from the databases
		std::string airfoil_path = "airfoil_databases/" + airfoil_type;
		std::ifstream airfoil_input(airfoil_path.c_str()); //in C++11, no need .c_str()
		
		if (airfoil_input.is_open()) {
			//------- READ AIRFOIL DATA AND APPEND TO ARRAY -------
			std::string airfoil_title, airfoil_coordinate;
			std::getline(airfoil_input, airfoil_title); // read the first line (title!)

			while (std::getline(airfoil_input, airfoil_coordinate)) {
				std::istringstream iss(airfoil_coordinate);
				double a, b;
				iss >> a >> b;
				x.push_back(a);
				y.push_back(b);
			}	

			max_node = x.size();
			
			// sharpen the airfoil's trailing edge, for kutta condition purpose.
			y[max_node - 1]	= 	0;
			y[0] 		= 	0;

			//reverse the index
			std::vector<double> x_temp	=	x;
			std::vector<double> y_temp	=	y;
			for (auto i = 1; i <= max_node; i++) {
				
				x[i - 1]	=	x_temp[max_node - i];
				y[i - 1]	=	y_temp[max_node - i];

			}

		} else {
			std::cout << "somehoew, the airfoil data is not exist!" << std::endl;
		}

	} else {
		std::cout << "Sorry, airfoil data is not found." << std::endl;
	}
}

bool Airfoil::airfoil_check_from_databases(std::string airfoil_input) {

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

std::vector<double> Airfoil::airfoil_element_length_calc(Airfoil_Parameters airfoil_pars, int max_node) {
	
	std::vector<double> length(max_node - 1); //the index is -1 from node index

	//make local variables
	std::vector<double> x	=	airfoil_pars.x;
	std::vector<double> y	=	airfoil_pars.y;

	for (auto i = 0; i < max_node - 1; i++) {

		length[i]	= 	math_f.length_two_points(x[i+1], x[i], y[i+1], y[i]);
	}
		
	return length;
}

std::vector<double> Airfoil::airfoil_beta_calc(Airfoil_Parameters airfoil_pars, int max_node) {

	std::vector<double> beta(max_node);

	//make local variables
	std::vector<double> x	=	airfoil_pars.x;
	std::vector<double> y	=	airfoil_pars.y;
	std::vector<double> s	=	airfoil_pars.s;

	double delta_x;
	double delta_y;
	double length;
	//edge first
	//calculate length of 2 neighborhood points
	length			= 	math_f.length_two_points(x[max_node-2], x[1], y[max_node-2], y[1]);
	beta[0]			=	2*M_PI - math_f.angle_cos_rule(length, s[max_node-2], s[0]);
	beta[max_node - 1]	=	beta[0];

	for (auto i = 1; i < max_node - 1; i++) {
		length		=	math_f.length_two_points(x[i-1], x[i+1], y[i-1], y[i+1]);
		beta[i]		=	2*M_PI - math_f.angle_cos_rule(length, s[i-1], s[i]);
	}
	
	return beta;
}

//calculate nx -> normal derivative on x direction
std::vector<double> Airfoil::airfoil_nx_calc(Airfoil_Parameters airfoil_pars, int max_node) {
	
	std::vector<double> nx(max_node);

	//make local variables
	std::vector<double> x	=	airfoil_pars.x;
	std::vector<double> y	=	airfoil_pars.y;
	std::vector<double> s	=	airfoil_pars.s;

	//first, calculate nx on the element (not on node!)
	//calculate delta x, delta_y, and temp_nx, temp_ny
	std::vector<double> delta_x(max_node - 1);
	std::vector<double> delta_y(max_node - 1);
	std::vector<double> temp_nx(max_node - 1);
	std::vector<double> temp_ny(max_node - 1);

	for (auto i = 0; i < max_node - 1; i++) {
		delta_x[i]	=	x[i+1] - x[i];
		delta_y[i]	=	y[i+1] - y[i];
		temp_nx[i]	=	-1*delta_y[i]/s[i];
		temp_ny[i]	=	delta_x[i]/s[i];
	}

	//we need to average the adjacent element into one node point.
	//for both edges (the same point obviously, just different node)
	nx[0]		=	math_f.avg_n(temp_nx[0], temp_nx[max_node-2], temp_ny[0], temp_ny[max_node-2]);
	nx[max_node-1]	=	nx[0];

	//for other nodes (i = 1 until i = max_node - 2)
	for (auto i = 1; i < max_node - 1; i++) {
		nx[i]		=	math_f.avg_n(temp_nx[i-1], temp_nx[i], temp_ny[i-1], temp_ny[i]);
	}

	return nx;
}

//calculate ny -> normal derivative on y direction
std::vector<double> Airfoil::airfoil_ny_calc(Airfoil_Parameters airfoil_pars, int max_node) {
	
	std::vector<double> ny(max_node);

	//make local variables
	std::vector<double> x	=	airfoil_pars.x;
	std::vector<double> y	=	airfoil_pars.y;
	std::vector<double> s	=	airfoil_pars.s;

	//first, calculate ny on the element (not on node!)
	//calculate delta x, delta_y, and temp_nx, temp_ny
	std::vector<double> delta_x(max_node - 1);
	std::vector<double> delta_y(max_node - 1);
	std::vector<double> temp_nx(max_node - 1);
	std::vector<double> temp_ny(max_node - 1);

	for (auto i = 0; i < max_node - 1; i++) {
		delta_x[i]	=	x[i+1] - x[i];
		delta_y[i]	=	y[i+1] - y[i];
		temp_nx[i]	=	-1*delta_y[i]/s[i];
		temp_ny[i]	=	delta_x[i]/s[i];
	}

	//we need to average the adjacent element into one node point.
	//for both edges (the same point obviously, just different node)
	ny[0]		=	math_f.avg_n(temp_ny[0], temp_ny[max_node-2], temp_nx[0], temp_nx[max_node-2]);
	ny[max_node-1]	=	ny[0];

	//for other nodes (i = 1 until i = max_node - 2)
	for (auto i = 1; i < max_node - 1; i++) {
		ny[i]		=	math_f.avg_n(temp_ny[i-1], temp_ny[i], temp_nx[i-1], temp_nx[i]);
	}

	return ny;
}

//function airfoil_eta_calc -> shape function calculations
std::vector<double> Airfoil::airfoil_eta_calc(Airfoil_Parameters airfoil_pars, int max_node) {

	std::vector<double> eta(max_node);

	//make local variables
	std::vector<double> s	=	airfoil_pars.s;

	//summation of the vector
	double sum_s	=	std::accumulate(s.begin(), s.end(), 0.0);

	//both edges
	eta[0]		=	-1;
	eta[max_node-1]	=	1;

	double temp;
	
	//for the rest of the nodes
	for (auto i = 1; i < max_node - 1; i++) {
		temp 	=	temp + s[i-1];
		eta[i]	=	(2*temp/sum_s) - 1;
	}

	return eta;
}

//calc N1
std::vector<double> Airfoil::airfoil_N1_calc(Airfoil_Parameters airfoil_pars, int max_node) {
	
	std::vector<double> N1(max_node);
		
	//make local vars
	std::vector<double> eta	=	airfoil_pars.eta;

	for (auto i = 0; i < max_node; i++) {
		N1[i]	=	0.5*(1 - eta[i]);
	}

	return N1;
}

//calc N1
std::vector<double> Airfoil::airfoil_N2_calc(Airfoil_Parameters airfoil_pars, int max_node) {
	
	std::vector<double> N2(max_node);
			
	//make local vars
	std::vector<double> eta	=	airfoil_pars.eta;

	for (auto i = 0; i < max_node; i++) {
		N2[i]	=	0.5*(1 + eta[i]);
	}
	
	return N2;
}
