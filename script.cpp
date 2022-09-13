#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "src/GMGSA.hpp"
#include "src/combinations.hpp"
#include "src/DRS.hpp"
#include "src/DPS.hpp"

#include "modelers/WigleyModeler.hpp"
#include "modelers/WigleyAnalyticMoments.hpp"

#include "util/WriteToFile.hpp"
#include "util/timer.hpp"

int main() {

	/*
	int N = 300;
	int ns = 1;
	int iter = 7;
	int omega = 1;
	timer t;
	t.begin();
	std::vector<std::vector<double>> T = DPS({ {-1,1},{-1,1} },
		std::vector<std::vector<double>>(),
		N,
		ns,
		iter,
		omega);
	t.display();
	std::cout << "s\n";
	
	std::string filename = std::to_string(N) + "_2D_DPS.dat";
	std::string message = std::to_string(N) + " samples, "
		+ std::to_string(ns) + " subpopulation size, "
		+ std::to_string(iter) + " iterations,  "
		+ std::to_string(omega) + "omega ";

	WriteToFile(T, filename, message);

	t.begin();
	std::vector<std::vector<double>> T_new = DPS({ {-1,1},{-1,1} },
		T,
		N,
		ns,
		iter,
		omega);
	t.display();
	std::cout << "s\n";

	filename = std::to_string(N) + "_next_2D_DPS.dat";
	message = "with input from previous sample generation, " + std::to_string(N) + " samples, "
		+ std::to_string(ns) + " subpopulation size, "
		+ std::to_string(iter) + " iterations,  "
		+ std::to_string(omega) + "omega ";

	WriteToFile(T_new, filename, message);

	system("gnuplot -p util/DPS_plot_temp.sh");
	
	return 0;
	*/

	/*
	double N = 324;
	double i = 1;

	std::vector<std::vector<double>> S = DRS({ {-1,1},{-1,1} }, std::vector<std::vector<double>>(), N, 1);

	//std::vector<std::vector<double>> S = DPS({ {-1,1},{-1,1} }, std::vector<std::vector<double>>(), N, 2, 3, 1);

	std::vector<std::vector<double>> S_new = DRS({ {-1,1},{-1,1} }, S, N, 7);

	//std::vector<std::vector<double>> S_new = DPS({ {-1,1},{-1,1} }, S, N, 2, 3, 1);
	
	WriteToFile(S, "temp1.dat");
	WriteToFile(S_new, "temp2.dat");
	system("gnuplot -p util/temp.sh");
	return 0;
	*/

	std::srand(6987879);

	WigleyModeler Wigley {{0.8,1.2},{0.08,0.12},{0.05,0.075},0.2,0,1};

	typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float< 30 > > BigFloat;
	WigleyAnalyticMoments<BigFloat> WigleyAnalytic{ Wigley };

	std::vector<std::vector<BigFloat>> T1, T2, T3;// Three parameters vs sample size
	std::vector<BigFloat> SI;//Sensitivity indices

	/*
	int order = 4;
	timer t;
	for (int N = 100; N <= 1000; N+=150) {
		t.begin();
		SI = GMGSA<WigleyAnalyticMoments<BigFloat>, BigFloat>(WigleyAnalytic, N, 4);
		T1.push_back({ double(N),SI.at(0),double(order) });
		T2.push_back({ double(N),SI.at(1),double(order) });
		T3.push_back({ double(N),SI.at(2),double(order) });

		std::cout << "Finished with N = " << N << " in ";
		t.display();
		std::cout << "s \n";
	}
	*/

	
	int N = 100;
	timer t;
	for (int order = 0; order <= 5; order += 1) {
		t.begin();
		SI = GMGSA<WigleyAnalyticMoments<BigFloat>, BigFloat>(WigleyAnalytic, N, order);
		T1.push_back({ double(order),SI.at(0),double(N) });
		T2.push_back({ double(order),SI.at(1),double(N) });
		T3.push_back({ double(order),SI.at(2),double(N) });

		std::cout << "Finished with order = " << order << " in ";
		t.display();
		std::cout << "s \n";
	}
	

	WriteToFile(T1, "Par1.dat");
	WriteToFile(T2, "Par2.dat");
	WriteToFile(T3, "Par3.dat");
	
	system("gnuplot -p util/script.sh");
	

	return 0;
}
