#include <iostream>
#include <vector>

#include "src/DPS.hpp"
#include "src/WriteToFile.hpp"

int main() {
	
	std::vector<std::vector<double>> X = {{0,1},{0,1},{0,1}};
	
	std::vector<std::vector<double>> S = DPS(X,std::vector<std::vector<double>>(),100,1,5,1);
	
	WriteToFile(S,"samples.dat");
	
	return 0;
}
