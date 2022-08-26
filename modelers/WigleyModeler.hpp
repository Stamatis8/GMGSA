#ifndef WIGLEYMODELER_HPP
#define WIGLEYMODELER_HPP

#include <vector>
#include <cmath>
#include <string>

#include "../src/smpl_triangulation/smpl_triangulation.hpp"

class WigleyModeler{
	/*
		This class was originally intended as a minimal template to be used in GMGSA.hpp. The model represents the Modified Wigley Hull-Form as
			in section 2 of:
			
			"A new mathematical hull‑form with 10‑shape parameters for evaluation of ship response in waves Sadaoki Matsui"
			
		Notes:
		
			- The domain has been expanded from [0,1]\times[0,1] to [-1,1]\times[-1,1] so as to include the entire hull
			
			- For vertical < 0 in the domain, the reflection of the point (x,y,z) = f(u,-v) is returned
			
			- For vertical == 0 in the domain, the centerline (ie y = 0, z = 0) is returned at the specified x
	*/

public:

	/* Constructor */
	
	WigleyModeler(
		std::vector<double> L_range,
		std::vector<double> B_range,
		std::vector<double> d_range,
		double c1_in = 0,
		double c2_in = 0,
		double c3_in = 0){
		
		this->c1 = c1_in;
		this->c2 = c2_in;
		this->c3 = c3_in;
		
		// L,B,D domains
		this->design_space_attribute = {L_range,B_range,d_range};
	
		// setting some initial design
		this->design = {
			(L_range.at(0)+L_range.at(1))/2,
			(B_range.at(0)+B_range.at(1))/2,
			(d_range.at(0)+d_range.at(1))/2};
	}

	/* Evaluate */
		
	std::vector<double> evaluate(std::vector<double> args){
		/*
			Description: Evaluates the modified wigley hull at xi,zeta for current design
			
			Inputs:
				- std::vector<double> args
					-1 <= args.at(0) <= 1
					-1 <= args.at(1) <= 1
			
			Output:
				- std::vector<double> point
					point on wigley hull according to conventions listed up top
					
		*/
		
		double xi = args.at(0);
		
		double zeta = args.at(1);
		
		std::vector<double> point(3,0);
		
		point.at(0) = this->design.at(0)*xi/2;
		
		double h = (1-std::pow(zeta,2))*(1-std::pow(xi,2))*(1 + this->c1*std::pow(xi,2) + this->c2*std::pow(xi,4))
				+ this->c3*std::pow(zeta,2)*(1-std::pow(zeta,8))*std::pow(1-std::pow(xi,2),4);
		
		if (zeta < 0){// reflect evaluate(xi,-zeta)
			zeta = -zeta;
			h = -h;
		}
		else if (zeta == 0){
			h = 0;
		}
		
		point.at(1) = this->design.at(1)*h/2;
		
		point.at(2) = zeta*this->design.at(2);
		
		return point;
	}

	/* Set Design */

	void set_design(std::vector<double> design_in){
		this->design = design_in;
	}
	
	/* Get Design Space */
	
	std::vector<std::vector<double>> design_space(){
		return this->design_space_attribute;
	}
	
	/* Get Domain */
	
	std::vector<std::vector<double>> domain(){
		return std::vector<std::vector<double>> {{-1,1},{-1,1}};
	}



protected:

	/* Attributes */

	std::vector<std::vector<double>> design_space_attribute;// ith parameter \in [design_space.at(i).at(0), design_space.at(i).at(1)]
	
	std::vector<double> design;// current_design
	
	double c1;// length
	
	double c2;// breadth
	
	double c3;// depth
	
};// WigleyModeler

#endif// WIGLEYMODELER_HPP
