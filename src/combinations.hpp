#ifndef COMBINATIONS_HPP
#define COMBINATIONS_HPP
	
#include <vector>
#include <cmath>

std::vector<std::vector<int>> combinations(int n, int m);

std::vector<std::vector<int>> combinations(int n, int m) {
	/*
		Description: Returns the m^n different ordered n-tuples of integers from 0 to m-1
			ie: n = 3, m = 2 -> m^n = 8 combinations
				( 0 , 0 , 0 )	( 1 , 0 , 0 )
				( 0 , 0 , 1 )	( 1 , 0 , 1 )
				( 0 , 1 , 0 )	( 1 , 1 , 0 )
				( 0 , 1 , 1 )	( 1 , 1 , 1 )
	*/

	auto pow = [](int base, int exponent) {// pow lambda
		int result = 1;
		for (int i = 0; i < exponent; i++) {
			result *= base;
		}
		return result;
	};

	int N = pow(m, n);// number of combinations

	// Initializing result to all zeros
	std::vector<std::vector<int>> c = std::vector<std::vector<int>>(N, std::vector<int>(n, 0));//combinations

	int row;// current row of c
	int dup, i_times;

	for (int col = n-1; col >= 0; col--) {// iterating through each column of comb, starting from the end

		dup = pow(m,n-col-1);// multiplicity of each element of {0,...,m-1}
		i_times = pow(m, col);// number of times to add {0,...,m-1}
		
		row = 0;
		for (int j = 0; j < i_times; j++) {// insert {0,...,m-1} std::pow(m,col) times with dup multiplicity
			for (int elem = 0; elem < m; elem++) {// insert {0,...,m-1} once with dup multiplicity
				for (int k = 0; k < dup; k++) {// insert elem dup times
					c.at(row).at(col) = elem;
					row++;
				}// j
			}// elem
		}//
	}// col

	return c;
}



#endif//COBINATIONS_HPP