/*
	Description: Perform GMGSA on a Wigley hull using analytic moment evaluation and the boost multiprecision library
		using N_samples for the sensitivity analysis and SSV of order SSV_order

		N_digits of precision are used
*/

#include <iostream>
#include <vector>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "../modelers/WigleyModeler.hpp"

#define GMGSA_USE_BOOST_MULTIPRECISION // must be defined before including WigleyAnalyticMoments.hpp and GMGSA.hpp
#include "../modelers/WigleyAnalyticMoments.hpp"

#include "../src/GMGSA.hpp"

int main() {

	WigleyModeler model{ {0.8,1.2},{0.08,0.12},{0.05,0.075},0.2,0,1 };// wigley model initialization

	const int N_digits = 30;// How many digits of precision to use

	typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float< N_digits > > 
		BigFloat;//defining arbitrary precision type
	
	WigleyAnalyticMoments<BigFloat> Wigley{ model };// attach analytic model calculator to wigley model

	std::vector<BigFloat> GSI; //generalized sensitivity indices

	int N_samples = 10;// number of samples
	int SSV_order = 9;// SSV order

	GSI = GMGSA<WigleyAnalyticMoments<BigFloat>,BigFloat>(Wigley, N_samples, SSV_order);

	for (int i = 0; i < GSI.size(); i++) {
		std::cout << "Sensitivity index for t" << i << " is: " << GSI.at(i) << std::endl;
	}

	return 0;
}