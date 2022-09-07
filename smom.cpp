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

	//WigleyAnalyticMoments<double> WigleyAnalytic{ Wigley };

	typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float< 80 > > MyFloat;
	WigleyAnalyticMoments<MyFloat> WigleyAnalytic{ Wigley };

	int p = 50;
	int q = 50;
	int r = 50;
	
	timer t;
	t.begin();
	std::cout << "Moment " << p << " - " << q << " - " << r << " is: " << WigleyAnalytic.moment(p, q, r) << std::endl;
	std::cout << "(calculated in ";
	t.display();
	std::cout << "s)" << std::endl;
	return 0;
}