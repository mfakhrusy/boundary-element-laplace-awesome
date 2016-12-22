#include <iostream>
#include <cmath>
#include <vector>

class Integral {

	public:
		double simpson(std::vector<double> x, std::vector<double> f_x);
};

class Test_Function {
	
	public:
		double func_1(double x);
};
int main() {

	Integral integral;
	Test_Function test_f;

	double x_initial	=	0;
	double x_final		=	10;
	int max_node		=	999999;
	std::vector<double> x(max_node);
	std::vector<double> f_x(max_node);
	for (auto i = 0; i < max_node; i++) {
		double zeta	=	static_cast<double>(i)/static_cast<double>(max_node-1);
		x[i]		=	x_initial*(1-zeta) + x_final*zeta;
		f_x[i]		=	test_f.func_1(x[i]);
	}

	//calculate integral
	double result		=	integral.simpson(x, f_x);
	std::cout << result << std::endl;
}

double Integral::simpson(std::vector<double> x, std::vector<double> f_x) {

	int max_node	=	x.size();
	double x_i	=	x[0];
	double x_f	=	x[max_node - 1];
	double temp	=	(x_f - x_i)/(3.0*static_cast<double>(max_node));

	double sum_even = 0;
	for (auto i = 2; i < max_node-1; i = i+2) {
		sum_even = sum_even + f_x[i];
	}

	double sum_odd = 0;
	for (auto i = 1; i < max_node-1; i = i+2) {
		sum_odd = sum_odd + f_x[i];
	}

	double result	=	temp*(f_x[0] + f_x[max_node-1] + 4*sum_odd + 2*sum_even);
	return result;
}

double Test_Function::func_1(double x) {	

	return pow(x,3);

}
