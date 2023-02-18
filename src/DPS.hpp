#ifndef DPS_HPP
#define DPS_HPP

#include <vector>
#include <cmath>
#include <algorithm>

#include "SimpleMonteCarlo.hpp"

std::vector<std::vector<double>> DPS(
	std::vector<std::vector<double>> X,
	std::vector<std::vector<double>> S_hat,
	int N,
	int ns,
	int i_max,
	double omega
	);
	
std::vector<std::vector<double>> DPS(
	std::vector<std::vector<double>> X,
	std::vector<std::vector<double>> S_hat,
	int N,
	int ns,
	int i_max,
	double omega
	)
{
	/*
		Dynamic Propagation Sampling (DPS)
			
		Source: Algorithm 2, from S. Khan, P. Kaklis - 'From regional sensitivity to intra-sensitivity for 
			parametric analysis of free-form shapes: Application to ship design'
		
		Input:
		
			- std::vector<std::vector<double>> X
				Design space limits. The ith parameter in X is bounded by [X.at(i-1).at(0), X.at(i-1).at(1)]
				X.size() >= 1
				X.at(i).size() == 2 for every i
			- std::vector<std::vector<double>> S_hat
				previous samples. Repulsion criterion, 'repulses' output from S_hat
			- int N
				Number of samples to construct
				N >= 1
			- int ns
				Number of designs in each sub-population p_L
				ns >= 1
			- int i_max
				Maximum number of iterations
				i_max >= 1
			- double omega
				Weight coefficient of Non-collapsing criterion
				
		Output:
		
			- std::vector<std::vector<double>> S
				 N samples from X, heuristically optimized to satisfy the three criteria listed in source
		
		Note: S_hat should be passed as an empty vector, if it is not to be used. Ie DPS(X,std::vector<std::vector<double>>(),N,ns,i_max,omega)
	*/

	bool UseSimpleMetaheuristic = true;
	//if true, attempts to generate new samples with simple monte carlo
	//if false, attempts to generate new samples according to derivative of criterion

	int n = X.size(); // number of parameters in design space
	std::vector<std::vector<std::vector<double>>> P(N, std::vector<std::vector<double>>(ns, std::vector<double>(n,0))); // population
	
	/* Creating Initial Population */
	
	for (int i = 0; i < N; i++){//sub-populations
		for (int j = 0; j < ns; j++){//designs
			for (int k = 0; k < n; k++){//design-parameters
				P.at(i).at(j).at(k) = X.at(k).at(0) + (X.at(k).at(1) - X.at(k).at(0))*rand()/RAND_MAX;
			}
		}
	}
	
	/* Initializing S by first entry of each sub-population	*/
	
	std::vector<std::vector<double>> S(N,std::vector<double>(n,0)); // final-samples
	
	for (int i = 0; i < N; i++){
		S.at(i) = P.at(i).at(0);
	}
	
	/* Creating Discretization of X (used in heuristic optimizer and F2 calculation) */ 
	
	std::vector<double> Delta(n,0);// discretization size of each parameter range (each is divided into N segments)
	for (int i = 0; i < n; i++){
		Delta.at(i) = (X.at(i).at(1)-X.at(i).at(0))/N;
	}
	std::vector<std::vector<double>> X_discrete(n,std::vector<double>(N+1,0));
	for (int i = 0; i < n; i++){
		for (int j = 0; j < (N+1); j++){
			X_discrete.at(i).at(j) = X.at(i).at(0) + (X.at(i).at(1)-X.at(i).at(0))*j/N;	
		}
	}
	
	/* Starting Iterations */
	
	std::vector<std::vector<double>> S_prime = S;
	
	double F1S, F2S, F3S;// S-criteria
	double F1S_prime, F2S_prime, F3S_prime;// S_prime-criteria 
	auto Enorm = [](std::vector<double> a, std::vector<double> b)// Euclidean norm lamda-expression
	{
		double norm = 0;
		for (int i = 0; i < a.size(); i++){
			norm += (a.at(i)-b.at(i))*(a.at(i)-b.at(i));
		}
		norm = std::sqrt(norm);
		return(norm);
	};
	double M;// Euclidean norm
	
	std::vector<double> x_c_prime(n,0);// heuristic x_c optimization
	std::vector<double> DF1(n,0);// derivative of F1S, used for meta-optimization
	std::vector<double> DF3(n,0);// derivative of F3S, used for meta-optimization
	std::vector<double> DF13(n,0);// derivative of F1S + F3S, used for meta-optimization
	
	for (int iteration = 0; iteration < i_max; iteration++){
		for (int L = 0; L < N; L++){
			std::vector<double> score(ns,0);// each element of P.at(L) will produce a score equal to F1+F2+F3
											// the element with the lowest score will be chosen for substitution
											// in S before moving to P.at(L+1)
			
			for (int c = 0; c < ns; c++){
			
				S.at(L) = P.at(L).at(c);
				
				/* Calculating F1S and DF1 */
				
				F1S = 0;
				DF1 = std::vector<double>(n,0);
				for (int p = 0; p < (N-1); p++){// Calculating F1S and DF1
					for (int q = p + 1; q < N; q++){
						M = Enorm(S.at(p),S.at(q));
						F1S += 1/(M*M);
						
						if(!UseSimpleMetaheuristic){
							if (p == L || q == L) {
								int coef = -2/(M*M*M*M);
								if (q == L) coef *= -1;
								for (int i = 0; i < n; i++){
									DF1.at(i) += coef*(S.at(p).at(i) - S.at(q).at(i));
								} 
							}
						}
					}
				}
				
				/* Calculating F3S and DF3 */
				
				if (S_hat.size() == 0) {// If no S_hat is passed, F3 criterion is ignored
					F3S = 0;
					DF3 =  std::vector<double>(n,0);
				}
				else {
					F3S = 0;
					DF3 = std::vector<double>(n,0);
					for (int p = 0; p < N; p++){
						for (int q = 0; q < S_hat.size(); q++){
							M = Enorm(S.at(p),S_hat.at(q));
							F3S += 1/(M*M);
							
							if(!UseSimpleMetaheuristic){
								if (p == L) {
									int coef = -2/(M*M*M*M);
									for (int i = 0; i < n; i++){
										DF3.at(i) += coef*(S.at(p).at(i) - S_hat.at(q).at(i));
									} 
								}
							}
						}
					}
				}
				
				/* Calculating meta-heuristically x_c_prime and updating S_prime */
				
				if(!UseSimpleMetaheuristic){
					for (int i = 0; i < n; i++){// Total effect of 1st and 3rd criterion
						DF13.at(i) = DF1.at(i) + DF3.at(i);
					}
					std::vector<double> DF13_abs = DF13;
					for (int i = 0; i < n; i++){
						if (DF13.at(i) < 0) DF13_abs.at(i) *= -1;
					}
					int i_maxchange = std::max_element(DF13_abs.begin(),DF13_abs.end()) - DF13_abs.begin();// parameter with max change on F1(S)
					
					x_c_prime = S.at(L);
					double delta = std::max(N*Delta.at(i_maxchange)/(5*(iteration+1)),Delta.at(i_maxchange)/2);// increment of i_maxchange parameter
																											// move less as iterations increase
					if (DF13.at(i_maxchange) < 0) {// increase i_maxchange parameter to decrease F13
						if (S.at(L).at(i_maxchange) + delta > X.at(i_maxchange).at(1)) {
							x_c_prime.at(i_maxchange) = X.at(i_maxchange).at(1); 
						}
						else {
							x_c_prime.at(i_maxchange) += delta;
						}
					}
					else {// decrease i_maxchange parameter to decrease F13
						delta *= -1;
						if (S.at(L).at(i_maxchange) + delta < X.at(i_maxchange).at(0)) {
							x_c_prime.at(i_maxchange) = X.at(i_maxchange).at(0); 
						}
						else {
							x_c_prime.at(i_maxchange) += delta;
						}
					}
				}
				else{
					x_c_prime = SimpleMonteCarlo(X,1).at(0);
				}
				
				S_prime = S;
				S_prime.at(L) = x_c_prime;
				
				/* Calculating F1S_prime */
				
				F1S_prime = 0;
				for (int p = 0; p < (N-1); p++){
					for (int q = p + 1; q < N; q++){
						M = Enorm(S_prime.at(p),S_prime.at(q));
						F1S_prime += 1/(M*M);
					}
				}
				
				/* Calculating F2S, F2S_prime */
				
				F2S = 0;
				F2S_prime = 0; 
				int i_closest = 0;// closest element of X_discrete to number under question
				std::vector<double> distances(N+1,0);// distances between number under question and elements of X_discrete
				for (int p = 0; p < (N-1); p++){
					for (int q = p+1; q < N; q++){
						for (int j = 0; j < n; j++){
						
							//F2 calculation for S
							if (std::abs(S.at(p).at(j) - S.at(q).at(j)) > Delta.at(j)){//Spj and Sqj cant be in the same discretized interval
																					   //if their distance is larger than said interval
								F2S += 0;
							}
							else {
								for (int i = 0; i < (N+1); i++){
									if ((S.at(p).at(j) - X_discrete.at(j).at(i)) < 0) {
										distances.at(i) = X_discrete.at(j).at(i) - S.at(p).at(j); 
									}
									else{
										distances.at(i) = S.at(p).at(j) - X_discrete.at(j).at(i);
									}
								}
								i_closest = std::min_element(distances.begin(),distances.end())-distances.begin();
								
								if ((S.at(p).at(j) - X_discrete.at(j).at(i_closest))*((S.at(q).at(j) - X_discrete.at(j).at(i_closest))) < 0){
									// if Spj, Sqj lie at oposite sides of node, they are at different intervals
									F2S += 0; 
								}
								else if (std::abs(S.at(q).at(j) - X_discrete.at(j).at(i_closest)) > Delta.at(j)){
									// if Sqj is further away than one interval from closest node to Spj and in the same side of it, they
									// cannot be at the same interval
									F2S += 0;
								}
								else{
									// Let node be the closest element of X_discrete.at(j) to Spj. Sqj, Spj are on the same interval iff 
									// Sqj is on the same side of node as is Spj and the distance between Sqj and node is less than one interval
									F2S += 1;
								}
							}
							
							//F2 calculation for S_prime
							if (std::abs(S_prime.at(p).at(j) - S_prime.at(q).at(j)) > Delta.at(j)){//Spj and Sqj cant be in the same discretized interval
																					   //if their distance is larger than said interval
								F2S_prime += 0;
							}
							else {
								for (int i = 0; i < (N+1); i++){
									if ((S_prime.at(p).at(j) - X_discrete.at(j).at(i)) < 0) {
										distances.at(i) = X_discrete.at(j).at(i) - S_prime.at(p).at(j); 
									}
									else{
										distances.at(i) = S_prime.at(p).at(j) - X_discrete.at(j).at(i);
									}
								}
								i_closest = std::min_element(distances.begin(),distances.end())-distances.begin();
								
								if ((S_prime.at(p).at(j) - X_discrete.at(j).at(i_closest))*((S_prime.at(q).at(j) - X_discrete.at(j).at(i_closest))) < 0){
									// if Spj, Sqj lie at oposite sides of node, they are at different intervals
									F2S_prime += 0; 
								}
								else if (std::abs(S_prime.at(q).at(j) - X_discrete.at(j).at(i_closest)) > Delta.at(j)){
									// if Sqj is further away than one interval from closest node to Spj and in the same side of it, they
									// cannot be at the same interval
									F2S_prime += 0;
								}
								else{
									// Let node be the closest element of X_discrete.at(j) to Spj. Sqj, Spj are on the same interval iff 
									// Sqj is on the same side of node as is Spj and the distance between Sqj and node is less than one interval
									F2S_prime += 1;
								}
							}
						}
					}
				}
				
				/* Calculating F3S_prime */
				
				if (S_hat.size() == 0) {// If no S_hat is passed, F3 criterion is ignored
					F3S_prime = 0;
				}
				else {
					F3S_prime = 0;
					for (int p = 0; p < N; p++){
						for (int q = 0; q < S_hat.size(); q++){
							M = Enorm(S_prime.at(p),S_hat.at(q));
							F3S_prime += 1/(M*M);
						}
					}
				}
				
				/* Comparing F, F_prime and updating P */
				score.at(c) = std::min(F1S + omega*F2S + F3S, F1S_prime + omega*F2S_prime + F3S_prime);		
				if (F1S + omega*F2S + F3S > F1S_prime + omega*F2S_prime + F3S_prime){
					P.at(L).at(c) = x_c_prime;
				}
					
			}//sub-population
			
			/* Updating Lth sample of S from best candidate of P.at(L) */
			
			S.at(L) = P.at(L).at(std::min_element(score.begin(),score.end()) - score.begin());
			
		}//population
		
	}//iteration
	
	return(S);
	
}// DPS function

#endif // DPS_HPP
