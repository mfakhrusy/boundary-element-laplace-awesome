#ifndef MAIN_LINEAR_SYSTEM_HPP
#define MAIN_LINEAR_SYSTEM_HPP

	class Main_linear_solver {
		private:
			//main 2 global constant
			std::vector<std::vector<double>> G_ij;
			std::vector<std::vector<double>> H_ij;

			//BC of main equation
			std::vector<double> G_iN_plus;
			std::vector<double> G_iN_minus;
			std::vector<double> del_phi_n;

			//BC also
			double phi_TE_tot;
			double del_phi_n_N_plus_inf;
			double del_phi_n_N_minus_inf;
			std::vector<double> phi_inf;

			//variables of each global constant (detailed)
			std::vector<std::vector<double>> g_ij_1;
			std::vector<std::vector<double>> g_ij_2;
			std::vector<std::vector<double>> h_ij_1;
			std::vector<std::vector<double>> h_ij_2;

			//BC of variables of the global constant
			std::vector<double> g_i_1;
			std::vector<double> h_i_1;

			std::vector<double> g_i_N;
			std::vector<double> g_i_N_plus_1;



		public:

	};
#endif
