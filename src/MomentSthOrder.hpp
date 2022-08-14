#ifndef MOMENTSTHORDER_HPP
#define MOMENTSTHORDER_HPP

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
	std::vector<double> a(3,0), b(3,0), c(3,0);// triangle is parametrized as (u,v) = au+bv+c
	double a_mag,b_mag;// Factors in analytical expression of I_t
	std::vector<double> n(3,0);// a,b exterior product
	int id; // Term in analytical expression of I_t
	
	
	for (int t = 0; t < triangles.size(); t++){
		
		/* Evaluating A through I */
		
		if (is_invariant){
		}
		else{
		}
		
		/* Evaluating I_t */
		
		
		M += I_t;
	}
	
	return(M);

}// MomentSthOrder()


#endif // MOMENTSTHORDER_HPP
