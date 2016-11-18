#ifndef READ_AIRFOIL_HPP
#define READ_AIRFOIL_HPP
	class Airfoil {
	
		private:
			bool check_airfoil_from_databases(std::string airfoil_input);

		public:
			void read_airfoil();
			int nodes;
			// x and y for airfoil coordinates, read it from upper trailing edge, rotate counter clockwise until lower trailing edge
			std::vector<double> x;
			std::vector<double> y;
	
	};
#endif
