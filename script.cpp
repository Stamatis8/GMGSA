#include <iostream>
#include <vector>

#include "src/GMGSA.hpp"

#include "modelers/WigleyModeler.hpp"
#include "modelers/WigleyAnalyticMoments.hpp"

#include "util/WriteToFile.hpp"
#include "util/timer.hpp"

int main() {
	
	WigleyModeler Wigley {{0.8,1.2},{0.08,0.12},{0.05,0.075},0.2,0,1};

	WigleyAnalyticMoments WigleyAnalytic {Wigley};
	
	std::vector<std::vector<double>> T1, T2, T3;// Three parameters vs sample size
	std::vector<double> SI;//Sensitivity indices

	int order = 4;
	
	timer t;
	for (int N = 100; N <= 500; N+=50) {
		t.begin();
		SI = GMGSA(WigleyAnalytic, N, 4);
		T1.push_back({ double(N),SI.at(0),double(order) });
		T2.push_back({ double(N),SI.at(1),double(order) });
		T3.push_back({ double(N),SI.at(2),double(order) });

		std::cout << "Finished with N = " << N << " in ";
		t.display();
		std::cout << "s \n";
	}

	WriteToFile(T1, "Par1.dat");
	WriteToFile(T2, "Par2.dat");
	WriteToFile(T3, "Par3.dat");
	
	system("gnuplot -p util/script.sh");

	return 0;
}
