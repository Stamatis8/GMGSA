#include <iostream>
#include <algorithm>
#include <vector>

#include "src/DPS.hpp"

int main() {
	
	std::vector<std::vector<double>> X = {{0,1},{0,1}};
	std::vector<std::vector<double>> S = DPS(X,std::vector<std::vector<double>>(),10,10,100,1);
	
	for (int i = 0; i<S.size();i++){
		for (int j = 0;j<S.at(i).size();j++){
			std::cout << S.at(i).at(j) << " ";
		}
		std::cout << std::endl;
	}
	return 1;
}
