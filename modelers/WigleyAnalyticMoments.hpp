#ifndef WIGLEYANALYTICMOMENTS_HPP
#define WIGLEYANALYTICMOMENTS_HPP

#include <vector>
#include <cmath>

#include "WigleyModeler.hpp"
#include "../src/geom_moments/NchooseK_cache.hpp"
template<typename scalar = double>
class WigleyAnalyticMoments: public WigleyModeler{

public:

	/* Constructor */

	WigleyAnalyticMoments(WigleyModeler m):WigleyModeler(m){}

	scalar moment(int p, int q, int r, bool is_translation_invariant = false, bool is_scaling_invariant = false){
		/*
			Description:
				  - Calculates the s = p + q + r order geometric moment of the current Wigley hull in modeler.
				  - is_translation_invariant  == true calculates the translation invariant of said moment
				  - is_scaling_invariant == true calculates scaling invariant of said moment
		*/
		
		if (p%2 != 0 || q%2 != 0){// due to longitudinal and transverse symmetry
			return 0;
		}
		
		scalar M;// moment
		scalar V;// volume

		scalar L = this->design.at(0);
		scalar B = this->design.at(1);
		scalar d = this->design.at(2);
		
		if (is_translation_invariant){// center design
			V = this->moment(0,0,0);
			scalar Cz = this->moment(0,0,1)/V;
			d = d - Cz;
			
			//x,y are already centered
		}				
		
		NchooseK_cache nk;// binomial coefficient cache
		
		scalar sum = 0;
		scalar Xi;
		scalar Zi;
		
		for (int i = 0; i <= q + 1; i++){
					
			/* Zi calculation */
			
			Zi = 0;
			for (int i1 = 0; i1 <= (q+1-i); i1++){
				for (int i2 = 0; i2 <= i; i2++){
					Zi += nk.get(q+1-i,i1)*nk.get(i,i2)
						*std::pow(-1,i1+i2)
						/(r+2*i+2*i1+8*i2+1);
				}// i2
			}// i1
			
			/* Xi calculation */

			Xi = 0;
			for (int i3 = 0; i3 <= (q+1+3*i); i3++){
				for (int i4 = 0; i4 <= (q+1-i); i4++){
					for (int i5 = 0; i5 <= i4; i5++){
						Xi += nk.get(q+1+3*i,i3)*nk.get(q+1-i,i4)*nk.get(i4,i5)
							*std::pow(-1,i3)
							*std::pow(this->c2,i4-i5)
							*std::pow(this->c1,i5)
							/(p+2*i3+4*i4-2*i5+1);
					}// i5
				}// i4
			}// i3
			
			/**/
		
			sum += nk.get(q+1,i)*std::pow(this->c3,i)*Xi*Zi;
		
		}// i
		
		M = 4
			*std::pow(L/2,p+1)
			*std::pow(B/2,q+1)
			*std::pow(d,r+1)
			*sum
			/(q+1);
			
		if (is_scaling_invariant){
			if (!is_translation_invariant){// calculate volume if it has not been calculated already
				V = this->moment(0,0,0);
			}
			
			M = M/(std::pow(V,(1+(p+q+r)/3)));
		}
		
		return M;	
		
	}
};

#endif //WIGLEYANALYTICMOMENTS_HPP
