#include <iostream>
#include <vector>

#include "src/GMGSA.hpp"
#include "src/SSV.hpp"
#include "modelers/WigleyModeler.hpp"
#include "modelers/stdModelerMoments.hpp"
#include "modelers/WigleyAnalyticMoments.hpp"

int main() {
	
	WigleyModeler Wigley {{80,120},{25,35},{17,19}};
	
	int order = 2;
	int order_analysis = 1;
	int N_samples = 10;
	
	stdModelerMoments<WigleyModeler> WigleyApprox {Wigley,5000};
	
	WigleyAnalyticMoments WigleyAna {Wigley};
	
	std::cout << WigleyAna.moment(2,2,0) << std::endl;
	std::cout << WigleyApprox.moment(2,2,0) << std::endl;
	return 0;
	//WigleyApprox.triangulate("Wigley.stl");
	
	std::vector<double> ssv = SSV(WigleyApprox,order);
	
	std::vector<double> SI = GMGSA(WigleyApprox,N_samples,order_analysis);
	
	std::vector<double> ssv_an = SSV(WigleyAna,order);
	
	std::vector<double> SI_an = GMGSA(WigleyAna,N_samples,order_analysis);
	
	std::cout << "SSV length difference: "<< std::endl;
	double norm=0;
	for(int i = 0; i < ssv_an.size();i++){
		norm += std::pow(std::abs(ssv_an.at(i)-ssv.at(i)),2);
	}
	norm = std::sqrt(norm);
	std::cout << norm << std::endl;
	return 0;
}
