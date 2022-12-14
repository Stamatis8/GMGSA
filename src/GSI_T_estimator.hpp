#ifndef GSI_T_ESTIMATOR_HPP
#define GSI_T_ESTIMATOR_HPP

#include <vector>
#include <cmath>

template<typename scalar>
scalar GSI_T_estimator(std::vector<std::vector<scalar>> Y, std::vector<std::vector<scalar>> Y_prime);

template<typename scalar>
scalar GSI_T_estimator(std::vector<std::vector<scalar>> Y, std::vector<std::vector<scalar>> Y_prime){
	/*
		Description: Estimates the Generalized Total Sensitivity Index (GSI_T) of some parameter unknown to the function. 
			The parameter information is contained in Y, Y_prime but is not needed. This is simply the implementation of the
			approximation closed form found in:
			
				"Sensitivity indices for multivariate outputs", Fabrice Gamboa, Alexandre Janon, Thierry Klein, Agn`es Lagnoux
				
			Which is a generalization to multivariate outputs of the second index found in
			
				"Asymptotic normality and efficiency of two Sobol index estimators", Alexandre Janon, Thierry Klein, 
					Agnes Lagnoux-Renaudie, Maëlle Nodet, Clémentine Prieur
					
		Input:
			std::vector<std::vector<double>> Y
				- Y.size() > 0
				- Y.size() must be equal to Y_prime.size()
				- Y.at(i).size() must be equal to Y.at(j).size() for every i,j
				- Y.at(i).size() must be equal to Y_prime.at(i).size() for every i
				
			std::vector<std::vector<double>> Y_prime
				- see Y
				
		Output:
			double 
				Sensitivity index approximation. See references
	*/
	
	scalar nominator = 0;
	scalar denominator = 0;
	
	scalar A, B, C; //nom = sum(A - B/N), denom = sum(C - B/N)
	
	int N = Y.size();
	int k = Y.at(0).size();
	
	for (int l = 0; l < k; l++){
		/* A,B,C */
		
		A = 0;
		B = 0;
		C = 0;
		for (int i = 0; i < N; i++){
			A += Y.at(i).at(l)*Y_prime.at(i).at(l);
			B += (Y.at(i).at(l) + Y_prime.at(i).at(l))/2;
			C += (Y.at(i).at(l) * Y.at(i).at(l) + Y_prime.at(i).at(l) * Y_prime.at(i).at(l))/2;
		}
		
		B = B*B;
		
		nominator += A - B/N;
		denominator += C - B/N;
	}
	
	return nominator/denominator;
}

#endif //GSI_T_ESTIMATOR_HPP
