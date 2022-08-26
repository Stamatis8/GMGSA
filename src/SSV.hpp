#ifndef SSV_HPP
#define SSV_HPP

#include <vector>

template<typename PM>
std::vector<double> SSV(PM modeler, int order);
	
template<typename PM>
std::vector<double> SSV(PM modeler, int order){
	/*
		Description: Constructs the Shape Signature Vector (SSV) of *order*, for *modeler*. SSV is
			defined in
			
			[1] "Geometric Moment-Dependent Global Sensitivity Analysis without Simulation Data: Application to Ship Hull Form Optimisation"
					Shahroz Khan, Panagiotis Kaklis, Andrea Serani, Matteo Diez, 2022
		
		PM typename requirements:
			PM.moment(int p, int q, int r, bool is_translation_invariant, bool is_scaling_invariant)
				- Calculates the s = p + q + r order geometric moment of the current design in modeler.
				- is_translation_invariant  == true calculates the translation invariant of said moment
				- is_scaling_invariant == true calculates scaling invariant of said moment
		
		Input:
			- PM modeler
				modeler class
			- int order
				order of Shape Signature Vector as in [1]
	*/
	
	std::vector<double> SSV;// Shape signature vector	
	int count;// counts elements in combinations
	
	for (int s = 0; s <= order; s++){
	
		if (s != 1){// s = 1 is ignored. See reference
	
			/* find all combinations i,j,k such that 0 <= i,j,k <= s and i+j+k = s*/
			
			std::vector<std::vector<int>> combinations((s+1)*(s+2)/2,std::vector<int>(3,0));
			count = 0;
			
			// i+j+k must equal s
			for (int i = 0; i <= s; i++){
				// j+k must equal s-i
				for (int j = 0; j <= s-i; j++){
					// k must equal s-i-j
					combinations.at(count) = {i, j, s-i-j};
					count++;
				}
			}
			
			/* evaluate each combination and add to SSV */
			for (int i = 0; i < combinations.size(); i++){
				if(s == 0){// save non-scaling invariant volume (which is equal to 1)
					SSV.push_back(modeler.moment(combinations.at(i).at(0),combinations.at(i).at(1),combinations.at(i).at(2),false,false));
				}
				else{
					SSV.push_back(modeler.moment(combinations.at(i).at(0),combinations.at(i).at(1),combinations.at(i).at(2),true,true));
				}
			}
		}
	}// s

	return SSV;

}// SSV()

#endif // SSV_HPP
