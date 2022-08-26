#ifndef PARAMETRICMODELER_HPP
#define PARAMETRICMODELER_HPP

#include <vector>

class ParametricModeler{
	/*
		This class is a minimal template. For more details see documentation
	*/

public:

	/* Constructor */
	
	ParametricModeler(){
		
	}
	
	/* Set Design */

	void set_design(std::vector<double> design_in){
		this->design = design_in;
	}
	
	/* Get Design Space */
	
	std::vector<std::vector<double>> design_space(){
		return this->design_space_attribute;
	}

	/* Evaluate */

	std::vector<double> evaluate(std::vector<double> args){
		
		std::vector<double> point;// evaluated model
		
		return point;
	}
	
	/* Get domain */
	
	std::vector<std::vector<double>> domain(){
		/* Evaluate and return domain for current design */
		
		std::vector<std::vector<double>> D;// domain for current design
		
		return D;
	}

private:

	/* Attributes */

	std::vector<std::vector<double>> design_space_attribute;// ith parameter \in [design_space.at(i).at(0), design_space.at(i).at(1)]

	std::vector<double> design;// current design
	
}; //ParametricModeler

#endif //PARAMETRICMODELER_HPP
