#ifndef SSV_HPP
#define SSV_HPP

#include <vector>

template<typename PM, typename scalar = double>
std::vector<scalar> SSV(PM modeler, int order);
	
template<typename PM, typename scalar = double>
std::vector<scalar> SSV(PM modeler, int order){
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

		Output:
			- std::vector<scalar> SSV
				contains all moments of order = 0, ..., 'order'
				behaviour can be changed based on defined macros (see below)

		Macros:
			- SSV_REMOVE_ZEROS
				- The moments in SSV that are *identically* zero are removed
			- SSV_EXACT_ORDER
				- Only moments whose order is exactly 'order' are included
	*/
	
	std::vector<scalar> SSV;// Shape signature vector	
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
			scalar M;//moment
			for (int i = 0; i < combinations.size(); i++){

#ifdef SSV_EXACT_ORDER// Do not calculate moments not equal to 'order'
				if (combinations.at(i).at(0) + combinations.at(i).at(1) + combinations.at(i).at(2) != order) continue;
#endif

				if(s == 0){// save non-scaling invariant volume (which is equal to 1)
					
#ifdef SSV_ALL_SCALING_INVARIANT
					M = modeler.moment(combinations.at(i).at(0), combinations.at(i).at(1), combinations.at(i).at(2), true, true);
#elif defined SSV_MOST_SCALING_INVARIANT// only non scaling invariant entry is volume when order == 0
					if (order == 0 || order == 1) {// if SSV is of order zero or one, calculate actual volume
						M = modeler.moment(combinations.at(i).at(0), combinations.at(i).at(1), combinations.at(i).at(2), false, false);
					}
					else {// if SSV is of order > 1, volume is set to 1 (scaling-invariant)
						M = modeler.moment(combinations.at(i).at(0), combinations.at(i).at(1), combinations.at(i).at(2), true, true);
					}
#else
					M = modeler.moment(combinations.at(i).at(0), combinations.at(i).at(1), combinations.at(i).at(2), false, false);
#endif

				}
				else{

#ifdef SSV_NONE_SCALING_INVARIANT
					M = modeler.moment(combinations.at(i).at(0), combinations.at(i).at(1), combinations.at(i).at(2), true, false);
#else
					M = modeler.moment(combinations.at(i).at(0), combinations.at(i).at(1), combinations.at(i).at(2), true, true);
#endif

				}

#ifdef SSV_REMOVE_ZEROS// Do not include moments identically equal to zero
				if (M == scalar(0)) continue;
#endif 

				SSV.push_back(M);
			}
		}
	}// s

	return SSV;

}// SSV()

#endif // SSV_HPP
