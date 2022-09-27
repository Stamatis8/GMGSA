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
			
			[1] "A new mathematical hull‑form with 10‑shape parameters for evaluation of ship response in waves Sadaoki Matsui"
			
		Notes:
		
			- The domain has been expanded from [0,1]\times[0,1] to [-1,1]\times[-1,1] so as to include the entire hull
			
			- For vertical < 0 in the domain, the reflection of the point (x,y,z) = f(u,-v) is returned
			
			- For vertical == 0 in the domain, the centerline (ie y = 0, z = 0) is returned at the specified x

		Parameters:

			Utilizing [1], there is more than one way to construct this modeler. The user will have the option to chose between
				these different parametrizations, based on the constructor overload they use. The constructor will always take 
				the following form:

				WigleyModeler(Length,Breadth,Depth,c1,c2,c3)

					** This is true in most cases. Ie see type = 2 below **

			Each of these arguments will be considered a parameter, if it is passed as a range of values instead of a singular one.

			In example:

				WigleyModeler MyWigley = WigleyModeler({10, 20},5,{6,7},0.1,0,0.3): 
					- Length is a parameter, whose range in the design-space is {10, 20}
					- Breadth is not a parameter, it is held fixed at 5
					- Depth is a parameter, whose range in the design-space is {6, 7}
					- c1 is not a parameter, it is held fixed at 0.1
					- c2 is not a parameter, it is held fixed at 0
					- c3 is not a parameter, it is held fixed at 0.3

			The various options are organized in the 'this->type' attribute. The options implemented so far are the following:

				type = 0:
					
					constructor: WigleyModeler({},{},{},c1,c2,c3)

					Parameters = Length, Breadth, Depth
					Non-Parameters = c1,c2,c3

				type = 1:

					constructor: WigleyModeler({},{},{},{},{},{})

					Parameters = Length, Breadth, Depth, c1 ,c2 ,c3

				type = 2:

					constructor: WigleyModeler(L,{},{},c1,c2,c3,"ratios")

					Parameters = B/L, d/L
					Non-Parameters = Length, c1, c2, c3

				type = 3:

					constructor: WigleyModeler(L,{},{},{},{},{},"ratios")

					Parameters = B/L, d/L, c1, c2, c3
					Non-Parameters = Length
				
				type = 4:

					constructor: WigleyModeler(L,B,d,{},{},{},"ratios")

					Parameters = c1, c2, c3
					Non-Parameters = Length, Breadth, Depth
			
	*/

public:

	/* Constructors */
	
	WigleyModeler(
		std::vector<double> L_range,
		std::vector<double> B_range,
		std::vector<double> d_range,
		double c1_in = 0,
		double c2_in = 0,
		double c3_in = 0){
		
		this->type = 0;// (see above)

		//Creating Design Space
		this->design_space_attribute = {L_range,B_range,d_range};
	
		//Initializing L,B,d,c1,c2,c3
		this->L = (L_range.at(0) + L_range.at(1)) / 2;
		this->B = (B_range.at(0) + B_range.at(1)) / 2;
		this->d = (d_range.at(0) + d_range.at(1)) / 2;
		this->c1 = c1_in; this->c2 = c2_in; this->c3 = c3_in;

		//Initializing design
		this->design = { this->L, this->B, this->d};
	}

	WigleyModeler(
		std::vector<double> L_range,
		std::vector<double> B_range,
		std::vector<double> d_range,
		std::vector<double> c1_range,
		std::vector<double> c2_range,
		std::vector<double> c3_range) {

		this->type = 1;// (see above)

		//Creating Design Space
		this->design_space_attribute = { L_range,B_range,d_range,c1_range,c2_range,c3_range };

		//Initializing L,B,d,c1,c2,c3
		this->L = (L_range.at(0) + L_range.at(1)) / 2;
		this->B = (B_range.at(0) + B_range.at(1)) / 2;
		this->d = (d_range.at(0) + d_range.at(1)) / 2;
		this->c1 = (c1_range.at(0) + c1_range.at(1)) / 2; 
		this->c2 = (c2_range.at(0) + c2_range.at(1)) / 2;
		this->c3 = (c3_range.at(0) + c3_range.at(1)) / 2;

		//Initializing design
		this->design = { this->L, this->B, this->d, this->c1, this->c2, this->c3 };
	}

	WigleyModeler(
		double L_in,
		std::vector<double> P1_range,
		std::vector<double> P2_range,
		double c1_in = 0,
		double c2_in = 0,
		double c3_in = 0,
		std::string r = "ratios") {

		if (r == "ratios") {

			this->type = 2;// (see above)

			//Creating Design Space
			this->design_space_attribute = { P1_range, P2_range };

			//Initializing L,B,d,c1,c2,c3
			this->L = L_in;
			double B_over_L = (P1_range.at(0) + P1_range.at(1)) / 2;
			double d_over_L = (P2_range.at(0) + P2_range.at(1)) / 2;
			this->B = B_over_L * this->L;
			this->d = d_over_L * this->L;
			this->c1 = c1_in; this->c2 = c2_in; this->c3 = c3_in;

			//Initializing design
			this->design = { B_over_L, d_over_L };
		}
	}

	WigleyModeler(
		double L_in,
		std::vector<double> P1_range,
		std::vector<double> P2_range,
		std::vector<double> c1_range,
		std::vector<double> c2_range,
		std::vector<double> c3_range,
		std::string r = "ratios") {

		if (r == "ratios") {

			this->type = 3;// (see above)

			//Creating Design Space
			this->design_space_attribute = { P1_range, P2_range, c1_range, c2_range, c3_range };

			//Initializing L,B,d,c1,c2,c3
			this->L = L_in;
			double B_over_L = (P1_range.at(0) + P1_range.at(1)) / 2;
			double d_over_L = (P2_range.at(0) + P2_range.at(1)) / 2;
			this->B = B_over_L * this->L;
			this->d = d_over_L * this->L;
			this->c1 = (c1_range.at(0) + c1_range.at(1)) / 2;
			this->c2 = (c2_range.at(0) + c2_range.at(1)) / 2;
			this->c3 = (c3_range.at(0) + c3_range.at(1)) / 2;

			//Initializing design
			this->design = { B_over_L, d_over_L, this->c1, this->c2, this->c3 };
		}
	}

	WigleyModeler(
		double L_in,
		double B_in,
		double d_in,
		std::vector<double> c1_range,
		std::vector<double> c2_range,
		std::vector<double> c3_range) {

		this->type = 4;// (see above)

		//Creating Design Space
		this->design_space_attribute = { c1_range, c2_range, c3_range };

		//Initializing L,B,d,c1,c2,c3
		this->L = L_in;
		this->B = B_in;
		this->d = d_in;
		this->c1 = (c1_range.at(0) + c1_range.at(1)) / 2;
		this->c2 = (c2_range.at(0) + c2_range.at(1)) / 2;
		this->c3 = (c3_range.at(0) + c3_range.at(1)) / 2;

		//Initializing design
		this->design = { this->c1, this->c2, this->c3 };
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
		
		double h = (1-std::pow(zeta,2))*(1-std::pow(xi,2))*(1 + this->c1*std::pow(xi,2) + this->c2*std::pow(xi,4))
				+ this->c3*std::pow(zeta,2)*(1-std::pow(zeta,8))*std::pow(1-std::pow(xi,2),4);
		
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

		switch (this->type) {
		case 0:// parameters are L,B,d

			this->L = this->design.at(0);
			this->B = this->design.at(1);
			this->d = this->design.at(2);

			break;
		case 1:// parameters are L,B,d,c1,c2,c3

			this->L = this->design.at(0);
			this->B = this->design.at(1);
			this->d = this->design.at(2);
			this->c1 = this->design.at(3);
			this->c2 = this->design.at(4);
			this->c3 = this->design.at(5);

			break;
		case 2:// parameters are B/L, d/L

			this->B = this->design.at(0) * this->L;
			this->d = this->design.at(1) * this->L;

			break;
		case 3:// parameters are B/L, d/L, c1, c2, c3

			this->B = this->design.at(0) * this->L;
			this->d = this->design.at(1) * this->L;
			this->c1 = this->design.at(2);
			this->c2 = this->design.at(3);
			this->c3 = this->design.at(4);

			break;

		case 4:

			this->c1 = this->design.at(0);
			this->c2 = this->design.at(1);
			this->c3 = this->design.at(2);

			break;
		}
	}
	
	/* Get Design Space */
	
	std::vector<std::vector<double>> design_space(){
		return this->design_space_attribute;
	}
	
	/* Get Domain */
	
	std::vector<std::vector<double>> domain(){
		return std::vector<std::vector<double>> {{-1,1},{-1,1}};
	}

	/* Type Attribute */

	unsigned int type;// type of modeler (see above)

	/* Get Particulars */

	std::vector<double> get_particulars() {
		return { this->L, this->B, this->d, this->c1, this->c2, this->c3 };
	}

protected:

	/* Attributes (private) */

	std::vector<std::vector<double>> design_space_attribute;// ith parameter \in [design_space.at(i).at(0), design_space.at(i).at(1)]
	
	std::vector<double> design;// current_design

	double L;// length

	double B;// breadth

	double d;// depth

	double c1;
	
	double c2;
	
	double c3;
	
};// WigleyModeler

#endif// WIGLEYMODELER_HPP
