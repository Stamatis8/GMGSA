#include <vector>
#include <iostream>

#include "../src/GMGSA.hpp"

#include "../modelers/WigleyModeler.hpp"
#include "../modelers/stdModelerMoments.hpp"

int main(){

	double c1 = 0.2;// length
	double c2 = 0;// breadth
	double c3 = 1;// depth
	int N_triangles = 200;// number of triangles to create mesh with

	WigleyModeler model { {0.8,1.2},{0.08,0.12},{0.05,0.075},0.2,0,1 };
	
	stdModelerMoments<WigleyModeler> Wigley { model, N_triangles };

	std::vector<double> GSI; //generalized sensitivity indices
	
	GSI = GMGSA(Wigley,50,2);// 50 samples, order 2 SSV
	
	for (int i = 0; i < GSI.size(); i++){	
		std::cout << "Sensitivity index for t" << i << " is: " << GSI.at(i) << std::endl;
	}

	return 0;
}
