/*
	Description: Perform GMGSA to a wigley hull using analytic moment calculation using N_samples for the sensitivity analysis and SSV of order
		SSV_order
*/
#include <iostream>
#include <vector>

#include "../modelers/WigleyModeler.hpp"
#include "../modelers/WigleyAnalyticMoments.hpp"

#include "../src/GMGSA.hpp"

int main() {

	WigleyModeler model{ {0.8,1.2},{0.08,0.12},{0.05,0.075},0.2,0,1 };// wigley model initialization

	WigleyAnalyticMoments<long double> Wigley { model };// attach analytic model calculator to wigley model

	std::vector<double> GSI; //generalized sensitivity indices

	int N_samples = 10000;// number of samples
	int SSV_order = 4;// SSV order

	GSI = GMGSA(Wigley, N_samples, SSV_order);

	for (int i = 0; i < GSI.size(); i++) {
		std::cout << "Sensitivity index for t" << i << " is: " << GSI.at(i) << std::endl;
	}

	return 0;
}