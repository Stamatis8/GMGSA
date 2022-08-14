#ifndef NCHOOSEK_HPP
#define NCHOOSEK_HPP

#include <vector>

namespace NchooseK {
	/*
		this namespace is used to access the precomputed nchoosek table throughout the program
	*/	
	
	struct table_struct{
		std::vector<std::vector<double>> nk;
		/*
			- the nk table contains (n,k) pairs. nk.at(i).at(j) is equal to (i choose j)
			- throughout the program execution, more and more values are calculated and included in nk
			- let M = nk.size. All elements nk.at(i).at(j), for i = 0,...,M-1 are calculated
		*/
		
		int get(int n,int k){
			/*
				Description: checks if (n,k) is calculated in nk. If it is calculated in nk, then nk.at(n).at(k) is returned.
					If it is not calculated then nk.size() - 1 <  n. Then, all values nk.at(i) for i = nk.size() - 1 to n are
					calculated
				
				Input:
					
					- int n
						
					- int k
						
				Output:
				
					- int out
						out = (n choose k)
						
				Note: Since at any poiny, all values nk.at(i) for i = 1,...,nk.size() are calculated, (n choose k) is calculated recursively
			*/
			
			/* Check if (n choose k) has been calculated */
			
			if ((nk.size()-1) >= n ) return(nk.at(n).at(k));
			
			/* If it has not been calculated, calculate it recursively */
			
			int prev = nk.size() - 1;// previously, maximum nk index was nk.at(prev)
		
			nk.resize(n+1);
			
			for (int i = prev+1; i <= n; i++){// nk.at(prev+1) through nk.at(n) must be filled
			
			}
		}
	}
}


#endif NCHOOSEK_HPP// NCHOOSEK_HPP
