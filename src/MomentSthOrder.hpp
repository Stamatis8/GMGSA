#ifndef MOMENTSTHORDER_HPP
#define MOMENTSTHORDER_HPP

#include <vector>
#include <cmath>
#include <algorithm>

#include "../util/NchooseK.hpp"

double MomentSthOrder(
	std::vector<std::vector<std::vector<double>>> triangles,
	int p,
	int q,
	int r,
	bool is_invariant
	);

double MomentSthOrder(
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
		
		n.at(0) = (a.at(1)*b.at(2) - a.at(2)*b.at(1))/(a_mag*b_mag); 
		n.at(1) = (a.at(2)*b.at(0) - a.at(0)*b.at(2))/(a_mag*b_mag);
		n.at(2) = (a.at(0)*b.at(1) - a.at(1)*b.at(0))/(a_mag*b_mag);
		
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
			volume = MomentSthOrder(triangles,0,0,0,false);
			
			/***** Todo: Handle volume == 0 *****/
			
			centroid.at(0) = MomentSthOrder(triangles,1,0,0,false)/volume;
			centroid.at(1) = MomentSthOrder(triangles,1,0,0,false)/volume;
			centroid.at(2) = MomentSthOrder(triangles,1,0,0,false)/volume;
			
			C -= centroid.at(0);
			F -= centroid.at(1);
			I -= centroid.at(2);
		}		
		
		/* Evaluating I_t */
		
		I_t = 0;
		
		for (int ip = 0; ip <= p; ip++){
			for (int iq = 0; iq <= q; iq++){
				for (int ir = 0; ir <= r; ir++){
				
					id = p + q + r + 1 - ip - iq - ir;
				
					/* Evaluating G1 G1_bar G1_tilde */
					G1 = 0;
					G1_bar = 0;
					G1_tilde = 0;
					
					for (int jp = 0; jp <= ip; jp++){
						for (int jq = 0; jq <= iq; jq++){
							for (int jr = 0; jr <= ir; jr++){
								for (int jd = 0; jd <= id + 1; jd++){
								
									power_row = std::pow(B,ip-jp)*std::pow(C,jp)
													  *std::pow(E,iq-jq)*std::pow(F,jq)
													  *std::pow(H,ir-jr)*std::pow(I,jr)*pow(-1,jd);//common factor in G1,G1_bar,G1_tilde
								
									if(jd != id + 1){// G1, G1_bar should be computed only in loop jd from 0 to id
										
										G1_bar += NK.get(ip,jp)*NK.get(iq,jq)*NK.get(ir,jr)*NK.get(id,jd)
											* power_row
											/ (p + q + r + 2 - jp - jq - jr -jd);
										
										G1_bar += NK.get(ip,jp)*NK.get(iq,jq)*NK.get(ir,jr)*NK.get(id,jd)
											* power_row
											/ ((p + q + r + 3 - jp - jq - jr -jd) * (p + q + r + 2 - jp - jq - jr -jd));
									}
									
									G1_tilde += NK.get(ip,jp)*NK.get(iq,jq)*NK.get(ir,jr)*NK.get(id+1,jd)
											* power_row
											/ (p + q + r + 3 - jp - jq - jr -jd);
								}//jd
							}//jr
						}//jq
					}//jp
				
					/*  Evaluating I_t */
					
					I_t += (n.at(0)*(B-A)/(p+1) + n.at(1)*(E-D)/(q+1) + n.at(2)*(H-G)/(r+1))*(G1-G1_bar)/id
						 + (n.at(0)*(C+A)/(p+1) + n.at(1)*(F+D)/(q+1) + n.at(2)*(I+G)/(r+1))*G1/id
						 - (n.at(0)*A/(p+1) + n.at(1)*D/(q+1) + n.at(2)*G/(r+1))*G1_tilde/(id*(id+1));
				
				}//ir
			}//iq
		}//ip
		I_t *= a_mag*b_mag/3;
		M += I_t;
	}//triangles
	
	if (is_invariant){//amending M (scaling invariance)
		M /= std::pow(volume,1+(p+q+r)/3);
	}
	
	return(M);

}// MomentSthOrder()


#endif // MOMENTSTHORDER_HPP
