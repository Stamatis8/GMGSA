#include <iostream>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cmath>
#include <vector>

#include "modelers/WigleyModeler.hpp"
#include "modelers/WigleyAnalyticMoments.hpp"
#include "util/WriteToFile.hpp"
#include "util/timer.hpp"

#include "src/pow_t.hpp"

int main() {

	WigleyModeler Wigley{ {0.8,1.2},{0.08,0.12},{0.05,0.075},0.2,0,1 };

	typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float< 80 > > MyFloat;
	WigleyAnalyticMoments<MyFloat> WigleyAnalytic{ Wigley };

	int p = 12;
	int q = 12;
	int r = 12;
	
	timer t;

	//Standard geometric moment
	t.begin();
	std::cout << "Moment " << p << " - " << q << " - " << r << " is: " << WigleyAnalytic.moment(p, q, r) << std::endl;
	std::cout << "(calculated in ";
	t.display();
	std::cout << "s)" << std::endl;

	//Translation invariant geometric moment
	t.begin();
	std::cout << "Translation invariant moment " << p << " - " << q << " - " << r << " is: " << WigleyAnalytic.moment(p, q, r, true, false) << std::endl;
	std::cout << "(calculated in ";
	t.display();
	std::cout << "s)" << std::endl;
	
	//Translation and scaling invariant geometric moment
	t.begin();
	std::cout << "Translation and scaling invariant moment " << p << " - " << q << " - " << r << " is: " << WigleyAnalytic.moment(p, q, r, true, true) << std::endl;
	std::cout << "(calculated in ";
	t.display();
	std::cout << "s)" << std::endl;

	return 0;
}