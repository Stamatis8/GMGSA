#include <iostream>
#include <vector>

#include "src/GMGSA.hpp"
#include "src/SSV.hpp"
#include "modelers/WigleyModeler.hpp"
#include "modelers/stdModelerMoments.hpp"

int main() {
	
	WigleyModeler Wigley {100,30,19};
	
	int order = 5;
	
	std::cout << "With approximate Wigley: " << std::endl;
	
	stdModelerMoments<WigleyModeler> WigleyApprox {Wigley,1800};
	
	WigleyApprox.set_design(std::vector<double> {0.5,0.3,0.2});
	WigleyApprox.triangulate("Wigley.stl");
	
	std::cout << "ssv:"<<std::endl;
	std::vector<double> ssv = SSV(WigleyApprox,order);
	for(int i = 0; i < ssv.size();i++){
		std::cout << ssv.at(i) << std::endl;
	}
	
	std::vector<double> SI = GMGSA(WigleyApprox,10,1);
	
	std::cout << "si:"<<std::endl;
	for(int i = 0; i < SI.size(); i++){
		std::cout << SI.at(i) << std::endl;
	}
	
	/*
	std::cout << "Volume is: " << WigleyApprox.moment(0,0,0) << std::endl;
	
	std::cout << "Centeroid is: ( "<< WigleyApprox.moment(1,0,0)/WigleyApprox.moment(0,0,0) 
		<< " , "<<
		WigleyApprox.moment(0,1,0)/WigleyApprox.moment(0,0,0) 
		<< " , "<<
		WigleyApprox.moment(0,0,1)/WigleyApprox.moment(0,0,0) 
		<< " )"<<std::endl;
	*/
	
	std::cout << "With analytic Wigley: " << std::endl;

	return 0;
}
