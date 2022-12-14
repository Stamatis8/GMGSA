#ifndef WIGLEYANALYTICMOMENTS_HPP
#define WIGLEYANALYTICMOMENTS_HPP

#include <vector>
#include <cmath>

#include "WigleyModeler.hpp"
#include "../src/geom_moments/NchooseK_cache.hpp"
#include "../src/pow_t.hpp"

#ifdef GMGSA_USE_BOOST_MULTIPRECISION// use the boost multiprecision library
									 // utilizes the boost::multiprecision::pow()
									 // function when exponent is not an integer
	#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

template<typename scalar>
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

		scalar L = this->L;
		scalar B = this->B;
		scalar d = this->d;
		
		scalar c = 0;// constant related to translation invariant formulation
		if (is_translation_invariant){
			 
			V = this->moment(0,0,0);
			c = this->moment(0, 0, 1) / (V * d);
			
			//x,y are already centered

		}				
		
		NchooseK_cache<scalar> nk;// binomial coefficient cache
		
		scalar sum = 0;
		scalar Xi;
		scalar Zi;

		for (int i = 0; i <= q + 1; i++){
					
			/* Zi calculation */
			
			Zi = 0;
			for (int i1 = 0; i1 <= (q+1-i); i1++){
				for (int i2 = 0; i2 <= i; i2++){
	
					if (is_translation_invariant) {
						// translation-invariant formulation

						for (int j = 0; j <= r; j++) {
							Zi = Zi + scalar(nk.get(q + 1 - i, i1))
								* scalar(nk.get(i, i2))
								* scalar(nk.get(r, j))
								* pow_t<scalar>(-1, i1 + i2 + j)
								* pow_t<scalar>(c,j)
								/ scalar(r - j + 2 * i + 2 * i1 + 8 * i2 + 1);
						}

					}
					else {// standard
						Zi = Zi + scalar(nk.get(q + 1 - i, i1))
							* scalar(nk.get(i, i2))
							* pow_t<scalar>(-1, i1 + i2)
							/ scalar(r + 2 * i + 2 * i1 + 8 * i2 + 1);
					}
					
				}// i2
			}// i1
			
			/* Xi calculation */

			Xi = 0;
			for (int i3 = 0; i3 <= (q+1+3*i); i3++){
				for (int i4 = 0; i4 <= (q+1-i); i4++){
					for (int i5 = 0; i5 <= i4; i5++){
						
						 Xi = Xi + scalar(nk.get(q + 1 + 3 * i, i3)) 
							* scalar(nk.get(q + 1 - i, i4))
							* scalar(nk.get(i4, i5))
							* pow_t<scalar>(-1, i3)
							* pow_t<scalar>(this->c2, i4 - i5)
							* pow_t<scalar>(this->c1, i5)
							/ scalar(p + 2 * i3 + 4 * i4 - 2 * i5 + 1);

					}// i5
				}// i4
			}// i3
			
			/**/
		
			sum = sum + scalar(nk.get(q + 1, i)) * pow_t<scalar>(this->c3, i) * Xi * Zi;

		}// i

		M = scalar(4)
			* pow_t<scalar>(L / 2, p + 1)
			* pow_t<scalar>(B / 2, q + 1)
			* pow_t<scalar>(d, r + 1)
			* sum
			/ scalar(q + 1);

		if (is_scaling_invariant){
			if (!is_translation_invariant){// calculate volume if it has not been calculated already
				V = this->moment(0,0,0);
			}

			#ifdef GMGSA_USE_BOOST_MULTIPRECISION
				M = M / (boost::multiprecision::pow(V, (1 + scalar(p + q + r) / 3)));
			#else
				M = M / (std::pow(V, (1 + scalar(p + q + r) / 3)));
			#endif
	
		}
		
		return M;	
		
	}
};

#endif //WIGLEYANALYTICMOMENTS_HPP
