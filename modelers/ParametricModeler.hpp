#ifndef PARAMETRICMODELER_HPP
#define PARAMETRICMODELER_HPP

#include <vector>

class ParametricModeler{
	/*
		This class is a minimal template to be used in GMGSA.hpp
	*/

public:

	/* Constructor */
	
	ParametricModeler(){
		
	}

	/* Triangulate */

	std::vector<std::vector<std::vector<double>>> triangulate(int N){
		/*
			Description: Construct a triangulation for this->design in approximately N triangles
		*/
		
		std::vector<std::vector<std::vector<double>>> triangulation;
		
		return triangulation;
	}

	/* Set Design */

	void set_design(std::vector<double> design_in){
		this->design = design_in;
	}
	
	/* Get Design Space */
	
	std::vector<std::vector<double>> design_space(){
		return this->design_space_attribute;
	}

private:

	/* Attributes */

	std::vector<std::vector<double>> design_space_attribute;// ith parameter \in [design_space.at(i).at(0), design_space.at(i).at(1)]

	std::vector<double> design;// current_design
	
}; //ParametricModeler

#endif //PARAMETRICMODELER_HPP
