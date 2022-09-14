/*
	Example: Perform GMGSA to a wigley hull modeler which approximates the geometric moments
		using a mesh of N_triangles using N_samples for the sensitivity analysis and SSV of order
		SSV_order
*/
#include <vector>
#include <iostream>

#include "../src/GMGSA.hpp"

#include "../modelers/WigleyModeler.hpp"
#include "../modelers/stdModelerMoments.hpp"

int main(){

	WigleyModeler model { {0.8,1.2},{0.08,0.12},{0.05,0.075},0.2,0,1 };// wigley model initialization
	
	int N_triangles = 200;// number of triangles to create mesh with

	stdModelerMoments<WigleyModeler> Wigley { model, N_triangles };// attach approximate moment calculator to model

	std::vector<double> GSI; //generalized sensitivity indices
	
	int N_samples = 50;// number of samples to use in sensitivity analysis
	int SSV_order = 2;// SSV order to use in sensitivity analysis 


	GSI = GMGSA(Wigley, N_samples, SSV_order);// Calculating GSI
	
	for (int i = 0; i < GSI.size(); i++){	
		std::cout << "Sensitivity index for t" << i << " is: " << GSI.at(i) << std::endl;
	}

	return 0;
}
