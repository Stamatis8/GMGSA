#ifndef GMGSA_HPP
#define GMGSA_HPP

#include <vector>

#include "GSI_T_estimator.hpp"
#include "DPS.hpp"
#include "DRS.hpp"
#include "SSV.hpp"
#include "geom_moments/geom_moments.hpp"

template<typename PM>
std::vector<double> GMGSA(PM modeler, int N, int order);

template<typename PM>
std::vector<double> GMGSA(PM modeler,int N, int order){
	/*
		Description: Performs Geometric Moment-dependent Global Sensitivity Analysis on modeler, according to
		
				[1] "Geometric Moment-Dependent Global Sensitivity Analysis without Simulation Data: Application to Ship Hull Form Optimisation"
					Shahroz Khan, Panagiotis Kaklis, Andrea Serani, Matteo Diez, 2022
			
			The output is the generalized total sensitivity indices
			
		class modeler requirements:
		
			- modeler.design_space()
				- returns an std::vector<std::vector<double>> vec
				- vec.at(i) is of size 2 and is equal to the range of the ith parameter
				- ie t_i \in [vec.at(i).at(0), vec.at(i).at(1)]
				- therefore it must also be true that number of parameters == vec.size()
				- modeler.design_space().size() > 0
			
			- modeler.set_design(std::vector<double> design)
				- design.size() == modeler.design_space().size()
				- design contains a set of parameters for the modeler
				- the modeler class must save this design until it is changed again
				
			( PM must be accepted by SSV() ):
			- modeler.moment(int p, int q, int r, bool is_translation_invariant, bool is_scaling_invariant)
				- Calculates the s = p + q + r order geometric moment of the current design in modeler.
				- is_translation_invariant  == true calculates the translation invariant of said moment
				- is_scaling_invariant == true calculates scaling invariant of said moment
		
		Input:
			- PM modeler
				- PM is a Parametric Modeler class, with the requirements listed above
			
			- int N
				- N > 0
				- number of samples to estimate the generalized sensitivity indices with	
			
			 - int order
				- order of shape signature vector
				
		Output:
			- std::vector<double> SI
				SI.at(i) is equal to the generalized total sensitivity index for the ith parameter
		
	*/
	
	bool use_DPS = false;// if true uses DPS.hpp for sample generation. Else uses DRS.hpp

	/* Generate N samples X */
	
	int sub_population_size = 2;// see references in DPS.hpp
	int max_iterations = 3;// see references in DPS.hpp
	double omega = 1;// see references in DPS.hpp

	std::vector<std::vector<double>> X;// samples

	if (use_DPS) {

		X = DPS(
			modeler.design_space(),
			std::vector<std::vector<double>>(),
			N,
			sub_population_size,
			max_iterations,
			omega
		);
	}
	else {
		X = DRS(
			modeler.design_space(),
			std::vector<std::vector<double>>(),
			N,
			1
		);
	}

	/* Using modeler and SSV calculate shape signature vector array Y */
	
	std::vector<std::vector<double>> Y(X.size(),std::vector<double>()); //Y.at(i) equals the SSV of X.at(i) of order
	std::vector<std::vector<std::vector<double>>> triangulation;// triangulation of design
	
	for (int design = 0; design < X.size(); design++){
		modeler.set_design(X.at(design));// choose design
		
		Y.at(design) = SSV(modeler,order);// calculate SSV for design
	}
	
	/* Initialize result SI and start calculation */
	
	std::vector<double> SI(modeler.design_space().size(),0);// SI.at(i) = ith sensitivity index

	/* Using DPS generate independent samples X_prime for ith parameter */

	std::vector<std::vector<double>> X_new;

	if (use_DPS) {

		X_new = DPS(
			modeler.design_space(),
			X,
			N,
			sub_population_size,
			max_iterations,
			omega
		);

	}
	else {

		X_new = DRS(
			modeler.design_space(),
			X,
			N,
			5
		);
	}

	std::vector<std::vector<double>> X_prime = X_new;
		// For each parameter, X_prime.at(i) will hold said parameter from X and set all others from X_new
	
	std::vector<std::vector<double>> Y_prime(X_prime.size(),std::vector<double>());
	
	for (int parameter = 0; parameter < modeler.design_space().size(); parameter++){// Calculate GSI total, for each parameter

		/* Calculating X_prime */
		
		for (int i = 0; i < N; i++){// only change parameter of ith sample in X_prime
			X_prime.at(i).at(parameter) = X.at(i).at(parameter);
		}

		/* Using modeler and SSV calculate Y_prime */
		
		for (int design = 0; design < X_prime.size(); design++){
			modeler.set_design(X_prime.at(design));// choose design
			
			Y_prime.at(design) = SSV(modeler,order);// calculate SSV for design
		}
		
		/* Using GSI_T_estimator, calculate SI.at(parameter) */
		
		SI.at(parameter) = GSI_T_estimator(Y,Y_prime);
		
		/* Undoing X_prime change to bring it back to X_prime = X_new */
		
		for (int i = 0; i < N; i++){// change parameter of ith sample in X_prime to X_new
			X_prime.at(i).at(parameter) = X_new.at(i).at(parameter);// Is it faster to simply assign X_prime = X_new each iteration?
		}
	}

	return SI;

}// GMGSA()


#endif //GMGSA_HPP
