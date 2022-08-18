#include <vector>
#include <iostream>

#include "../src/GMGSA.hpp"

#include "../modelers/WigleyModeler.hpp"

int main(){

	double L = 100;// length
	double B = 30;// breadth
	double d = 15;// depth

	WigleyModeler model(L,B,d);
	
	std::vector<double> GSI; //generalized sensitivity indices
	
	GSI = GMGSA(model,50,200,2);// 50 samples, 200 triangles per sample, order 2 SSV
	
	for (int i = 0; i < GSI.size(); i++){	
		std::cout << "Sensitivity index for t" << i << " is: " << GSI.at(i) << std::endl;
	}

	return 0;
}
