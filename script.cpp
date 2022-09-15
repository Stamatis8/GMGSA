#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>

#include <boost/multiprecision/cpp_dec_float.hpp>

//#define GMGSA_USE_BOOST_MULTIPRECISION
#define GMGSA_USE_TBB

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

	std::srand(12322);

	std::vector<double> L_range = { 0.8,1.2 };
	std::vector<double> B_range = { 0.08,0.12 };
	std::vector<double> T_range = { 0.05,0.075 };
	double c1 = 0.2;
	double c2 = 0;
	double c3 = 1;

	WigleyModeler Wigley { L_range,B_range,T_range,c1,c2,c3};

	//typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float< 30 > > BigFloat;

	typedef long double BigFloat;

	WigleyAnalyticMoments<BigFloat> WigleyAnalytic{ Wigley };

	std::vector<std::vector<BigFloat>> T1, T2, T3;// Three parameters vs sample size
	std::vector<BigFloat> SI;//Sensitivity indices

	std::string title, xlabel, ylabel;
	title = "Wigley Hull with L,B,T parameters\\n"
		" L in [" + std::to_string(L_range.at(0)) + "," + std::to_string(L_range.at(1)) + "], B in [" + std::to_string(B_range.at(0)) + "," + std::to_string(B_range.at(1)) + "], T in [" + std::to_string(T_range.at(0)) + "," + std::to_string(T_range.at(1)) + "]\\n"
		" c1 = " + std::to_string(c1) + ", c2 =" + std::to_string(c2) + ", c3 = " + std::to_string(c3) + "\\n";

	/*
	int order = 6;
	timer t;
	for (int N = 1000; N <= 5000; N+=1000) {
		t.begin();
		SI = GMGSA<WigleyAnalyticMoments<BigFloat>, BigFloat>(WigleyAnalytic, N, order);
		T1.push_back({ double(N),SI.at(0),double(order) });
		T2.push_back({ double(N),SI.at(1),double(order) });
		T3.push_back({ double(N),SI.at(2),double(order) });

		std::cout << "Finished with N = " << N << " in ";
		t.display();
		std::cout << "s" << std::endl;
	}
	*/

	
	int N = 10000;
	timer t;
	for (int order = 0; order <= 6; order += 1) {
		t.begin();
		
		SI = GMGSA<WigleyAnalyticMoments<BigFloat>, BigFloat>(WigleyAnalytic, N, order);
	
		T1.push_back({ double(order),SI.at(0),double(N) });
		T2.push_back({ double(order),SI.at(1),double(N) });
		T3.push_back({ double(order),SI.at(2),double(N) });

		std::cout << "Finished with order = " << order << " in ";
		t.display();
		std::cout << "s " << std::endl;
	}
	title += "using simple MonteCarlo sampling with " + std::to_string(N) + " samples\\n";
	xlabel = "order";
	ylabel = "SI";
	

	WriteToFile(T1, "Par1.dat");
	WriteToFile(T2, "Par2.dat");
	WriteToFile(T3, "Par3.dat");
	
	//system("gnuplot -p util/script.sh");

	std::string message = 
		"set title \"" + title + "\"; set xlabel \"" + xlabel + "\"; set ylabel \"" + ylabel + "\" rotate by 0;"
		"set key outside; set key right;"
		"plot 'Par1.dat' with linespoint title \"L\", 'Par2.dat' with linespoint title \"B\", 'Par3.dat' with linespoint title \"T\";"
		;
	
	WriteToFile(message,"gnuplot.sh");
	
	return 0;
}
