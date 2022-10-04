#include <vector>
#include <cmath>
#include <string>

#include "WigleyModeler.hpp"
#include "../src/pow_t.hpp"
#include "../src/geom_moments/NchooseK_cache.hpp"
#include "../src/combinations.hpp"

template<typename scalar>
std::vector<std::vector<scalar>> WigleyAnalyticSobol(WigleyModeler modeler, double order) {
/*
	Description:

		[1] introduces a generalization to multivariate outputs of the second sensitivity index investigated
		in [2], the so-called Sobol sensitivity index. With the help of [3], this function calculates
		analytically said sensitivity indexes for the SSV (see [4]) of the wigley hull family in section 2 of [5]

		By default, SI from all orders up to 'order' are calculated and for all parameters

	Input:

		WigleyModeler modeler
			- wigley modeler to perform sensitivity analysis on 
			- To determine which are the parameters of interest (ie L,B,T or B/L, T/L or etc...) 
				the .type attribute of WigleyModeler is utilized

				**IMPORTANT** only L,B,T is treated as of now

		double order
			- order of SSV to calculate

	Output:

		std::vector<std::vector<scalar>> out
			- out.at(i) = Sensitivity indices for order i+1 if i>0
			- out.at(i).at(j) = jth parameter sensitivity index for order i+1 if i>0

	References:

		[1] Gamboa, Sensitivity Indices for multivariate outputs, 2013

		[2] A. Janon, Asymptotic normality and efficiency of two sobol index estimators, 2012

		[3] I. M. Sobol, Sensitivity estimates for nonlinear mathematical models Math. 
			Modeling Comput Experiment, 1(4):407�414 (1995), 1993.

		[4] S. Khan, P Kaklis, A Serani, M Diez, \emph{Geometric Moment-Dependent Global Sensitivity Analysis without  
			Simulation Data: Application to Ship Hull Form Optimisation}, Computer-Aided Design 151 (2022) 103339

		[5] Sadaoki Matsui (2022), A new mathematical hull-form with 10-shape parameters for 
			evaluation of ship response in waves, Journal of Marine Science and Technology 27:508�521

*/

	/* Calculation of number of elements in SSV */

	double s_bar = 1;//number of elements in SSV

	for (int i = 2; i <= order; i++) {
		s_bar += (i + 1)*(i + 2) / 2;
	}

	/* Calculation of (p_j,q_j,r_j) for each index j in SSV */

	std::vector<std::vector<double>> combinations(s_bar, std::vector<double>(3, 0));
		//combinations.at(j) = p_j,q_j,r_j

	combinations.at(0) = { 0,0,0 };//Initializing s = 0

	double c_index = 1;// combinations index

	for (double s = 2; s <= order; s++) {
		// i+j+k must equal s
		for (double i = 0; i <= s; i++) {
			// j+k must equal s-i
			for (double j = 0; j <= s - i; j++) {
				// k must equal s-i-j
				combinations.at(c_index) = { i, j, s - i - j };
				c_index++;
			}
		}
	}

	WriteToFile(combinations, "DEBUG_combinations");//DEBUG

	/* Retrieving modeler constant parameters and design space */

	std::vector<double> particulars = modeler.get_particulars();// L,B,T,c1,c2,c3 of modeler
	//scalar L = particulars.at(0);
	//scalar B = particulars.at(1);
	//scalar d = particulars.at(2);
	scalar c1 = particulars.at(3);
	scalar c2 = particulars.at(4);
	scalar c3 = particulars.at(5);
	
	std::vector<std::vector<double>> D = modeler.design_space();// design space
	
	scalar L_0 = D.at(0).at(0);	scalar L_1 = D.at(0).at(1);
	scalar B_0 = D.at(1).at(0);	scalar B_1 = D.at(1).at(1);
	scalar d_0 = D.at(2).at(0);	scalar d_1 = D.at(2).at(1);

	scalar DeltaL = L_1 - L_0;
	scalar DeltaB = B_1 - B_0;
	scalar Deltad = d_1 - d_0;

	scalar c = 3 * (scalar(1 / 3) + c1 / 15 + c2 / 35 + c3 * 128 / 945) / (scalar(1 / 3) + c1 / 15 + c2 / 35 + c3 * 256 / 3465) / 4;
		//c is an fj related constant (see formulation)

	/* Initialization of loop through combinations for j \in [2,s_bar] */ 
	
	std::vector<std::vector<scalar>> out = std::vector<std::vector<scalar>>(order+1,std::vector<scalar>(3,0));// output (see above)
	double out_index = 2;// out.at(out_index) will be filled next
	double j_next = 6;// at j=j_next, S^u as far as it is calculated will be added to out.at(out_index)

	double p_j, q_j, r_j;//MI^{p_j,q_j,r_j} is jth element of SSV
	double P_j, Q_j, R_j;//see formulation

	scalar fj;// fj = f(c1,c2,c3) at p_j,q_j,r_j
	scalar Zi, Xi, sum;// fj related

	scalar IjL, IjB, Ijd;// see formulation
	scalar Ij2L, Ij2B, Ij2d;// see formulation

	scalar cj, NjL, NjB, Njd, Dj;// see formulation

	scalar nominatorL = 0,nominatorB = 0,nominatord = 0,denominator = 0;// 

	NchooseK_cache<scalar> nk;// binomial coefficient cache

	for (double j = 2; j < s_bar; j++) {

		p_j = combinations.at(j).at(0);
		q_j = combinations.at(j).at(1);
		r_j = combinations.at(j).at(2);

		P_j = 2 * p_j - q_j - r_j;
		Q_j = 2 * q_j - p_j - r_j;
		R_j = 2 * r_j - p_j - q_j;

		/* Calculate fj */

		sum = 0;

		for (int i = 0; i <= q_j + 1; i++) {

			/* Zi calculation */

			Zi = 0;
			for (int i1 = 0; i1 <= (q_j + 1 - i); i1++) {
				for (int i2 = 0; i2 <= i; i2++) {
					for (int j = 0; j <= r_j; j++) {
						Zi = Zi + scalar(nk.get(q_j + 1 - i, i1))
							* scalar(nk.get(i, i2))
							* scalar(nk.get(r_j, j))
							* pow_t<scalar>(-1, i1 + i2 + j)
							* pow_t<scalar>(c, j)
							/ scalar(r_j - j + 2 * i + 2 * i1 + 8 * i2 + 1);
					}
				}// i2
			}// i1

			/* Xi calculation */

			Xi = 0;
			for (int i3 = 0; i3 <= (q_j + 1 + 3 * i); i3++) {
				for (int i4 = 0; i4 <= (q_j + 1 - i); i4++) {
					for (int i5 = 0; i5 <= i4; i5++) {

						Xi = Xi + scalar(nk.get(q_j + 1 + 3 * i, i3))
							* scalar(nk.get(q_j + 1 - i, i4))
							* scalar(nk.get(i4, i5))
							* pow_t<scalar>(-1, i3)
							* pow_t<scalar>(c2, i4 - i5)
							* pow_t<scalar>(c1, i5)
							/ scalar(p_j + 2 * i3 + 4 * i4 - 2 * i5 + 1);

					}// i5
				}// i4
			}// i3

			/**/

			sum = sum + scalar(nk.get(q_j + 1, i)) * pow_t<scalar>(c3, i) * Xi * Zi;

		}// i

		fj = sum /
			(std::pow(4 * (scalar(1 / 3) + c1 / 15 + c2 / 35 + c3 * 256 / 3465) / 3, (p_j + q_j + r_j) / 3 + 1) *
				std::pow(2, p_j + q_j) * (q_j + 1));

		/* Calculate IjL,IjB,Ijd,Ij2L,Ij2B,Ij2d */

		if (P_j == -3) {
			IjL = std::log(L_1 / L_0);
		}
		else {
			IjL = (std::pow(L_1, P_j / 3 + 1) - std::pow(L_0, P_j / 3 + 1)) * 3 / (P_j + 3);
		}

		if (Q_j == -3) {
			IjB = std::log(B_1 / B_0);
		}
		else {
			IjB = (std::pow(B_1, Q_j / 3 + 1) - std::pow(B_0, Q_j / 3 + 1)) * 3 / (Q_j + 3);
		}

		if (R_j == -3) {
			Ijd = std::log(d_1 / d_0);
		}
		else {
			Ijd = (std::pow(d_1, R_j / 3 + 1) - std::pow(d_0, R_j / 3 + 1)) * 3 / (R_j + 3);
		}

		if (2*P_j == -3) {
			Ij2L = std::log(L_1 / L_0);
		}
		else {
			Ij2L = (std::pow(L_1, P_j * 2 / 3 + 1) - std::pow(L_0, P_j * 2 / 3 + 1)) * 3 / (P_j * 2 + 3);
		}

		if (2*Q_j == -3) {
			Ij2B = std::log(B_1 / B_0);
		}
		else {
			Ij2B = (std::pow(B_1, Q_j * 2 / 3 + 1) - std::pow(B_0, Q_j * 2 / 3 + 1)) * 3 / (Q_j * 2 + 3);
		}

		if (2*R_j == -3) {
			Ij2d = std::log(d_1 / d_0);
		}
		else {
			Ij2d = (std::pow(d_1, R_j * 2/ 3 + 1) - std::pow(d_0, R_j * 2 / 3 + 1)) * 3 / (R_j * 2 + 3);
		}

		/* Calculate cj */

		cj = fj * IjL * IjB * Ijd / (DeltaL * DeltaB * Deltad);
		/* Calculate Dj */

		Dj = fj*fj*Ij2L*Ij2B*Ij2d/(DeltaL*DeltaB*Deltad) - cj*cj;

		/* Calculate Nj for all parameters */

		NjL = Ij2L*std::pow(fj*IjB*Ijd/(DeltaB*Deltad),2)/DeltaL - cj*cj;

		NjB = Ij2B * std::pow(fj * IjL * Ijd / (DeltaL * Deltad), 2) / DeltaB - cj * cj;

		Njd = Ij2d * std::pow(fj * IjL * IjB / (DeltaL * DeltaB), 2) / Deltad - cj * cj;

		/* Add Nj, Dj to global Nominator, Denominator */

		denominator += Dj;
		nominatorL += NjL;
		nominatorB += NjB;
		nominatord += Njd;

		/* Save index to out and continue */
		
		if(j==j_next){
			out.at(out_index) = {nominatorL/denominator,nominatorB/denominator,nominatord/denominator};
			out_index++;
			j_next += (out_index+1)*(out_index+2)/2;
		}	
	}

	return out;
}