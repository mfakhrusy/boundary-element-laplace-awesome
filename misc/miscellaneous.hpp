#ifndef MISCELLANEOUS_HPP
#define MISCELLANEOUS_HPP

class Miscellaneous {
	
	void print_to_file(std::vector<double> f, std::string filename);
	void print_to_file(std::vector<std::vector<double>> f, std::string filename);

	void neumann_bc(std::vector<double> &F);

	//friend classes below can access the functions above (which are privates)
	friend class Matrix_Solver;
	friend class Matrix_Init;
	friend class Post_Computation;
};

#endif
