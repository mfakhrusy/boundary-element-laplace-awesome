#include "../global.hpp"
#include "miscellaneous.hpp"

void Miscellaneous::print_to_file(std::vector<double> f, std::string filename) {

	std::ofstream output_file;
	output_file.open("output-temp/" + filename);
	for (auto i = 0; i < f.size(); i++) {
		output_file << i << " " << f[i];
		output_file << std::endl;
	}
	output_file.close();

}

void Miscellaneous::print_to_file(std::vector<std::vector<double>> f, std::string filename) {

	std::ofstream output_file;
	output_file.open("output-temp/" + filename);
	for (auto i = 0; i < f.size(); i++) {
		for (auto j = 0; j < f[i].size(); j++) {
			output_file << f[i][j] << " ";
		}
		output_file << std::endl;
	}
	output_file.close();

}
