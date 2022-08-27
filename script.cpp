#include <iostream>
#include <vector>

#include "src/GMGSA.hpp"
#include "modelers/WigleyModeler.hpp"
#include "modelers/WigleyAnalyticMoments.hpp"

int main() {
	
	WigleyModeler Wigley {{0.8,1.2},{0.08,0.12},{0.05,0.075},0.2,0,1};

	WigleyAnalyticMoments WigleyAnalytic {Wigley};
	
	std::vector<double> SI = GMGSA(WigleyAnalytic,200,4);
	
	for (int i = 0; i < SI.size(); i++){
		std::cout << SI.at(i) << std::endl;
	}

	return 0;
}
