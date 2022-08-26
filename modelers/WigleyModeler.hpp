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
	
	WigleyModeler(double L_in, double B_in, double d_in){
		this->L = L_in;
		this->B = B_in;
		this->d = d_in;
		
		// c1,c2,c3 domains
		this->design_space_attribute = {{0,1},{0,1},{0,1}};
	
		// original wigley hull
		this->design = {0,0,0};
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
		
		point.at(0) = this->L*xi/2;
		
		double h = (1-std::pow(zeta,2))*(1-std::pow(xi,2))*(1 + this->design.at(0)*std::pow(xi,2) + this->design.at(1)*std::pow(xi,4))
				+ this->design.at(2)*std::pow(zeta,2)*(1-std::pow(zeta,8))*std::pow(1-std::pow(xi,2),4);
		
		if (zeta < 0){// reflect evaluate(xi,-zeta)
			zeta = -zeta;
			h = -h;
		}
		else if (zeta == 0){
			h = 0;
		}
		
		point.at(1) = this->B*h/2;
		
		point.at(2) = zeta*this->d;
		
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

private:

	/* Attributes */

	std::vector<std::vector<double>> design_space_attribute;// ith parameter \in [design_space.at(i).at(0), design_space.at(i).at(1)]

	std::vector<double> design;// current_design
	
	double L;// length
	
	double B;// breadth
	
	double d;// depth
	
};// WigleyModeler

#endif// WIGLEYMODELER_HPP
