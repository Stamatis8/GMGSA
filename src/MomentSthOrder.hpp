#ifndef MOMENTSTHORDER_HPP
#define MOMENTSTHORDER_HPP

#include <vector>
#include <cmath>
#include <algorithm>

#include "../util/NchooseK.hpp"

double MomentSthOrder(
	std::vector<std::vector<std::vector<double>>> triangles,
	int i,
	int j,
	int k,
	int r,
	bool is_invariant
	);

double MomentSthOrder(
	std::vector<std::vector<std::vector<double>>> triangles,
	int i,
	int j,
	int k,
	int r,
	bool is_invariant
	)
{
	/*
		Description: The input variable triangles describes the triangulation of a surface S, which encloses a volume G. This function calculates the
			s = (p+q+r) order geometric moment of G, according to Pozo's et.al. - *Approximate Series Algorithm For Geometric Moments* which can be found in
			section 4 of [1]. Further, the translation and scaling invariant moments of G can also be calculated as in section 3.2.1 of [2] 
					
		References:
		
			[1] Efficient 3D Geometric and Zernike moments computation from unstructured surface meshes J. M. Pozo, M. C. Villa-Uriol,
				A. F. Frangi, Senior Member, IEEE
				
			[2] Geometric Moment-Dependent Global Sensitivity Analysis without Simulation Data: Application to Ship Hull Form Optimisation, S. Khan, 
				P. Kaklis, 2022
		
		Input:
		
			- std::vector<std::vector<std::vector<double>>> triangles
				triangles.at(p) is the pth triangle with vertices triangles.at(p).at(q), q=0,1,2. Each vertex must have three elements
			- int i
				i >= 0
			- int j
				j >= 0
			- int k
				k >= 0
			- int r
				0 <= r <= i+j+k
			- bool is_translation_invariant
				if true, calculates the translation invariant moment form of G (ie centered at it's center of volume)
				By definition all first order moments are equal to zero, so if i+j+k = 1, set M = 0
			- bool is_scaling_invariant
				if true, calculates the scaling invariant moment form of G 
				By definition volume is equal to 1, so if i+j+k = 0, set M = 1
		Output:
		
			- double M
				the s = (i+j+k) order moment of G, according to options is_translation_invariant and is_scaling_invariant
	*/

	/* If is_translation_invariant, translate mesh */
	
	/* Initializing variables */
	
	int T = triangles.size();// number of triangles
	
	std::vector<std::vector<double>> p_bar = // p_bar.at(i) = centroid of ith triangle
		std::vector<std::vector<double>>(T,std::vector<double>(3,0)); 
	
	std::vector<std::vector<double>> u_hat = // Sec.4 of [2]
		std::vector<std::vector<double>>(T,std::vector<double>(3,0));
	
	std::vector<std::vector<double>> v_hat = // Sec.4 of [2]
		std::vector<std::vector<double>>(T,std::vector<double>(3,0));
	
	//area c lamda expression
	
	double lamda = 0;
	
	
	/* If is_scaling_invariant, divide final moment M by volume^(1+s/3) */

	return 0;
}

double DEBUG_TEMP(
		double A, 
		double B, 
		double C, 
		double D,
		double E,
		double F,
		double G,
		double H,
		double I,
		int a,
		int b,
		int c){
		double result = 0;
		int id = 0;
		
		NchooseK_cache NK;// N choose K cache
		NK.get(std::max(a,std::max(b,c)),0);// Precomputing N choose k pairs for all k and N <= max(p,q,r)
		
		for (int ia = 0; ia <= a; ia++){
			for (int ib = 0; ib <= b; ib++){
				for (int ic = 0; ic <= c; ic++){
				
					id = a + b + c + 1 - ia - ib - ic;
				
					for (int ja = 0; ja <= ia; ja++){
						for (int jb = 0; jb <= ib; jb++){
							for (int jc = 0; jc <= ic; jc++){
								for (int jd = 0; jd <= id; jd++){
									result += NK.get(a,ia)*NK.get(b,ib)*NK.get(c,ic)*NK.get(ia,ja)*NK.get(ib,jb)*NK.get(ic,jc)*NK.get(id,jd)
											* std::pow(A,a-ia)*std::pow(D,b-ib)*std::pow(G,c-ic)
											* std::pow(B,ia-ja)*std::pow(C,ja)*std::pow(E,ib-jb)
											* std::pow(F,jb)*std::pow(H,ic-jc)*std::pow(I,jc)
											* std::pow(-1,jd)/(id*(a+b+c+2-ja-jb-jc-jd));
								}//jd
							}//jc
						}//jb
					}//ja
					
				}//ic
			}//ib
		}// ia
		
		return result;
}// DEBUG_TEMP()


double old_MomentSthOrder(
	std::vector<std::vector<std::vector<double>>> triangles,
	int p,
	int q,
	int r,
	bool is_invariant
	)
{
	/*
		Description: The input variable triangles describes the triangulation of some surface S, which encloses some volume G. This function calculates the
			s = (p+q+r) order moment of G, as described in the appendix of 
			
			[1] Geometric Moment-Dependent Global Sensitivity Analysis without Simulation Data: Application to Ship Hull Form Optimisation - S. Khan - P. Kaklis -
			2022
		
		Input:
		
			- std::vector<std::vector<std::vector<double>>> triangles
				triangles.at(i) is the ith triangle with vertices triangles.at(i).at(j), j=0,1,2. Each vertex must have three elements
			- int p
				p >= 0
			- int q
				q >= 0
			- int r
				r >= 0
			- bool is_invariant
				calculates the invariant moment form of G (invariant to scaling and translation transformations) as described in eq.8 of [1] 
	
		Output:
		
			- double M
				the s = (p+q+r) moment of G
				
		Note: The method used in this implementation is the analytical evaluation of the integral of f (see appendix of [1]) over each triangle.
	*/
	
	double M = 0;
	double I_t = 0;// The integral of f (see appendix of [1]) over the triangle t
	
	double A,B,C,D,E,F,G,H,I;// Factors in analytitcal expression of I_t
	double G1,G1_bar,G1_tilde;// Summation terms in I_t expression
	double power_row;//common factor in G1,G1_bar,G1_tilde	
	
	std::vector<double> a(3,0), b(3,0), c(3,0);// triangle is parametrized as (u,v) = au+bv+c
	double a_mag,b_mag;// Factors in analytical expression of I_t
	std::vector<double> n(3,0);// a,b exterior product
	int id; // Term in analytical expression of I_t
	
	NchooseK_cache NK;// N choose K cache
	NK.get(std::max(p,std::max(q,r)),0);// Precomputing N choose k pairs for all k and N <= max(p,q,r)
	
	double volume = 0;// In case invariant M is requested
	std::vector<double> centroid(3,0);// In case invariant M is requested
	
	for (int t = 0; t < triangles.size(); t++){
		
		/* Evaluating A through I, a, b, c, n */
		
		a.at(0) = triangles.at(t).at(1).at(0) - triangles.at(t).at(0).at(0);
		a.at(1) = triangles.at(t).at(1).at(1) - triangles.at(t).at(0).at(1);
		a.at(2) = triangles.at(t).at(1).at(2) - triangles.at(t).at(0).at(2);
		
		b.at(0) = triangles.at(t).at(2).at(0) - triangles.at(t).at(0).at(0);
		b.at(1) = triangles.at(t).at(2).at(1) - triangles.at(t).at(0).at(1);
		b.at(2) = triangles.at(t).at(2).at(2) - triangles.at(t).at(0).at(2);
		
		c = triangles.at(t).at(0); 
		
		a_mag = std::sqrt(a.at(0)*a.at(0) + a.at(1)*a.at(1) + a.at(2)*a.at(2));
		b_mag = std::sqrt(b.at(0)*b.at(0) + b.at(1)*b.at(1) + b.at(2)*b.at(2));
		
		//n.at(0) = (a.at(1)*b.at(2) - a.at(2)*b.at(1))/(a_mag*b_mag); 
		//n.at(1) = (a.at(2)*b.at(0) - a.at(0)*b.at(2))/(a_mag*b_mag);
		//n.at(2) = (a.at(0)*b.at(1) - a.at(1)*b.at(0))/(a_mag*b_mag);
		
		n.at(0) = (a.at(1)*b.at(2) - a.at(2)*b.at(1)); 
		n.at(1) = (a.at(2)*b.at(0) - a.at(0)*b.at(2));
		n.at(2) = (a.at(0)*b.at(1) - a.at(1)*b.at(0));
		double n_mag = std::sqrt(n.at(0)*n.at(0) + n.at(1)*n.at(1) + n.at(2)*n.at(2));
		n.at(0) /= n_mag;
		n.at(1) /= n_mag;
		n.at(2) /= n_mag;

		
		A = a.at(0);
		B = b.at(0);
		C = c.at(0);
		D = a.at(1);
		E = b.at(1);
		F = c.at(1);
		G = a.at(2);
		H = b.at(2);
		I = c.at(2);
	
		if (is_invariant){//amending A through I (translational invariance)
			volume = old_MomentSthOrder(triangles,0,0,0,false);
			
			/***** Todo: Handle volume == 0 *****/
			
			centroid.at(0) = old_MomentSthOrder(triangles,1,0,0,false)/volume;
			centroid.at(1) = old_MomentSthOrder(triangles,1,0,0,false)/volume;
			centroid.at(2) = old_MomentSthOrder(triangles,1,0,0,false)/volume;
			
			C -= centroid.at(0);
			F -= centroid.at(1);
			I -= centroid.at(2);
		}		
		
		/* Evaluating I_t */
		
		/* DEBUG: Evaluating I_t another way */
		
		double term1 = 0;
		double term2 = 0;
		double term3 = 0;
		
		term1 = DEBUG_TEMP(A,B,C,D,E,F,G,H,I,p+1,q,r);
		term2 = DEBUG_TEMP(A,B,C,D,E,F,G,H,I,p,q+1,r);
		term3 = DEBUG_TEMP(A,B,C,D,E,F,G,H,I,p,q,r+1);
		
		term1 *= n.at(0)/(p+1);
		term2 *= n.at(1)/(q+1);
		term3 *= n.at(2)/(r+1);
		
		//I_t = a_mag*b_mag*(term1 + term2 + term3)/3;
		
		double ab = a.at(0)*b.at(0) + a.at(1)*b.at(1) + a.at(2)*b.at(2);
		I_t = std::sqrt(a_mag*a_mag*b_mag*b_mag - ab*ab)*(term1+term2+term3)/3;
				
		M += I_t;
	}//triangles
	
	if (is_invariant){//amending M (scaling invariance)
		M /= std::pow(volume,1+(p+q+r)/3);
	}
	
	return M;

}// MomentSthOrder()

#endif // MOMENTSTHORDER_HPP
