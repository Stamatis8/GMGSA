#ifndef SSV_HPP
#define SSV_HPP

#include <vector>
#include <string>

#include "MomentSthOrder.hpp"
#include "stl2vec.hpp"

std::vector<double> SSV(
	std::vector<std::vector<std::vector<double>>> triangles, int order);

std::vector<double> SSV(std::string filename, int order);

std::vector<double> SSV(
	std::vector<std::vector<std::vector<double>>> triangles, int order){
	/*
		Description: Constructs the Shape Signature Vector (SSV) of order of the volume enclosed by triangles. SSV is
			defined in
			
			"Geometric Moment-Dependent Global Sensitivity Analysis without Simulation Data: Application to Ship Hull Form Optimisation"
					Shahroz Khan, Panagiotis Kaklis, Andrea Serani, Matteo Diez, 2022
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
					SSV.push_back(MomentSthOrder(triangles,combinations.at(i).at(0),combinations.at(i).at(1),combinations.at(i).at(2),s,false,false));
				}
				else{
					SSV.push_back(MomentSthOrder(triangles,combinations.at(i).at(0),combinations.at(i).at(1),combinations.at(i).at(2),s,true,true));
				}
			}
		}
	}// s

	return SSV;

}// SSV()


std::vector<double> SSV(std::string filename, int order){
	/*
		Description: Extracts triangles from filename and passes them to SSV()
	*/
	
	return SSV(stl2vec(filename),order);
}
#endif // SSV_HPP
