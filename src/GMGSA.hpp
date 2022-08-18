#ifndef GMGSA_HPP
#define GMGSA_HPP

#include <vector>

#include "GSI_T_estimator.hpp"
#include "DPS.hpp"
#include "geom_moments/geom_moments.hpp"

template<class PM>
std::vector<double> GMGSA(PM modeler, int N ,int N_triangles, int order);

template<class PM>
std::vector<double> GMGSA(PM modeler,int N, int N_triangles, int order){
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
				
			- modeler.triangulate(N_triangles)
				- returns a triangulation of the currently set design in *approximately* N_triangles
				- the triangulation is returned in std::vector<std::vector<std::vector<double>>> triang format
				- the triangulation must be oriented
				- triang.at(i) is the ith triangle
				- triang.at(i).at(j) is the jth vertex of the ith triangle ( 0 <= j <= 2 )
				- triang.at(i).at(j).at(k) is the kth coordinate of the jth vertex of the ith triangle ( 0 <= k <= 2 )
		
		Input:
			- PM modeler
				- PM is a Parametric Modeler class, with the requirements listed above
			
			- int N
				- N > 0
				- number of samples to estimate the generalized sensitivity indices with	
			
			 - int N_triangles
				- N_triangles > 2
				- *indicator* for number of triangles to estimate the shape signature vectors with
			
			 - int order
				- order of shape signature vector
				
		Output:
			- std::vector<double> SI
				SI.at(i) is equal to the generalized total sensitivity index for the ith parameter
		
	*/
	
	/* Generate N samples X using DPS */
	
	int sub_population_size = 2;// see references in DPS.hpp
	int max_iterations = 5;// see references in DPS.hpp
	double omega = 1;// see references in DPS.hpp
	
	std::vector<std::vector<double>> X = DPS(
		modeler.design_space(),
		std::vector<std::vector<double>>(),
		N,
		sub_population_size,
		max_iterations,
		omega
	);
	
	/* Using modeler and SSV calculate shape signature vector array Y */
	
	std::vector<std::vector<double>> Y(X.size(),std::vector<double>()); //Y.at(i) equals the SSV of X.at(i) of order
	std::vector<std::vector<std::vector<double>>> triangulation;// triangulation of design
	
	for (int design = 0; design < X.size(); design++){
		modeler.set_design(X.at(design));// choose design
		
		triangulation = modeler.triangulate(N_triangles);// triangulate design
		
		Y.at(design) = SSV(triangulation,order);// calculate SSV for design
	}
	
	/* Initialize result SI and start calculation */
	
	std::vector<double> SI(modeler.design_space().size(),0);// SI.at(i) = ith sensitivity index
	
	/* Using DPS generate independent samples X_prime for ith parameter */
		
	std::vector<std::vector<double>> X_new = DPS(// Generate X_new with repulsion criterion (see references in DPS.hpp) applied to X
		modeler.design_space(),
		X,
		N,
		sub_population_size,
		max_iterations,
		omega
	);
	
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
		
			triangulation = modeler.triangulate(N_triangles);// triangulate design
		
			Y_prime.at(design) = SSV(triangulation,order);// calculate SSV for design
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
