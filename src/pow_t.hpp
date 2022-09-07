#ifndef POW_T_HPP
#define POW_T_HPP

template<typename scalar>
scalar pow_t(scalar base, int exponent) {
	/*
		Description: returns base^exponent. As with std::pow(), 0^0 is defined as 1

		typename scalar requirements:
			- scalar(0) constructor
			- scalar(1) constructor
			- * operator which returns scalar (multiplication)
			- + operator which returns scalar (addition)
	*/

	if (base == scalar(0)){
		if (exponent != 0) {
			return 0; // 0^a, a > 0 is 0
		}
		else {
			return 1; // std::pow() convention 0^0 is 1
		}
	}
	else if (exponent > 0) {
		return base * pow_t<scalar>(base, exponent - 1);
	}
	else {// exponent == 0 and base != 0
		return scalar(1);
	}
}

#endif// POW_T_HPP