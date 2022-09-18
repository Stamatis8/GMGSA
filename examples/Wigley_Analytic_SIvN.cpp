/*
	Description: Perform Sensitivity analysis for the Wigley hull for a fixed SSV order and N in a certain range. After how many N, do
		the sensitivity indices converge? Capability for multiple runs with different random seeding
*/

#include <vector>
#include <iostream>
#include <string>

//Float Multiprecision Library
//#define GMGSA_USE_BOOST_MULTIPRECISION
#ifdef GMGSA_USE_BOOST_MULTIPRECISION
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

//Thread Building Blocks
#define GMGSA_USE_TBB

//SSV Options
//#define SSV_REMOVE_ZEROS
//#define SSV_EXACT_ORDER

#include "../src/GMGSA.hpp"
#include "../modelers/WigleyModeler.hpp"
#include "../modelers/WigleyAnalyticMoments.hpp"

#include "../util/timer.hpp"
#include "../util/WriteToFile.hpp"

int main() {

	/* Experiment initialization */

	std::vector<long int> N_range = { 2000, 20000 };// SSV order range of experiment
	long int N_step = 2000;// calculate at N_range.at(0) + i*N_step until N_range.at(1)
	int order = 0;// number of samples to be used in experiment
	int runs = 10;// number of runs of experiment

	std::srand(2321);// random seed

	/* Wigley modeler initialization */

	std::vector<double> L_range = { 0.8,1.2 };
	std::vector<double> B_range = { 0.08,0.12 };
	std::vector<double> T_range = { 0.05,0.075 };
	double c1 = 0.2;
	double c2 = 0;
	double c3 = 1;

	WigleyModeler model{ L_range, B_range, T_range, c1, c2, c3 };

	/* Wigley moment calculator initialization */

#ifdef GMGSA_USE_BOOST_MULTIPRECISION
	typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float< 30 > > BigFloat;
#else
	typedef long double BigFloat;
#endif

	WigleyAnalyticMoments<BigFloat> Wigley{ model };

	/* Calculation */

	std::vector<std::vector<BigFloat>> T1, T2, T3;// Three parameter data
	std::vector<BigFloat> SI;//Sensitivity indices

	std::string title, xlabel, ylabel;
	title = "Wigley Hull with L,B,T parameters\\n"
		" L in [" + std::to_string(L_range.at(0)) + "," + std::to_string(L_range.at(1)) + "], B in [" + std::to_string(B_range.at(0)) + "," + std::to_string(B_range.at(1)) + "], T in [" + std::to_string(T_range.at(0)) + "," + std::to_string(T_range.at(1)) + "]\\n"
		" c1 = " + std::to_string(c1) + ", c2 =" + std::to_string(c2) + ", c3 = " + std::to_string(c3) + "\\n";

	timer t;
	for (int run = 1; run <= runs; run++) {
		std::cout << "Starting run (" << std::to_string(run) << " / " << std::to_string(runs) << ")" << std::endl;

		for (int N = N_range.at(0); N <= N_range.at(1); N += N_step) {
			t.begin();

			SI = GMGSA<WigleyAnalyticMoments<BigFloat>, BigFloat>(Wigley, N, order);

			T1.push_back({ double(N),SI.at(0),double(order) });
			T2.push_back({ double(N),SI.at(1),double(order) });
			T3.push_back({ double(N),SI.at(2),double(order) });

			std::cout << "Finished with N = " << N << " in ";
			t.display();
			std::cout << "s " << std::endl;
		}

		WriteToFile(T1, "data/Par1_" + std::to_string(run) + ".dat");
		WriteToFile(T2, "data/Par2_" + std::to_string(run) + ".dat");
		WriteToFile(T3, "data/Par3_" + std::to_string(run) + ".dat");

		T1.clear();T2.clear();T3.clear();

		std::cout << "Finished run (" << std::to_string(run) << " / " << std::to_string(runs) << ")" << std::endl << std::endl;
	}


	/* Constructing gnuplot script */

	title += "using simple MonteCarlo sampling. SSV order =  " + std::to_string(order) + " \\n";
	xlabel = "N";
	ylabel = "SI";

	std::string message =
		"set title \"" + title + "\"; set xlabel \"" + xlabel + "\"; set ylabel \"" + ylabel + "\" rotate by 0;"
		"set key outside;"
		;

	std::string filename1, filename2, filename3;
	for (int run = 1; run <= runs; run++) {
		filename1 = "data/Par1_" + std::to_string(run) + ".dat";
		filename2 = "data/Par2_" + std::to_string(run) + ".dat";
		filename3 = "data/Par3_" + std::to_string(run) + ".dat";

		if (run == 1) {
			message += "plot '" + filename1 + "' with linespoint title \"L\" lc rgb \"dark-violet\" pt 1, '" + filename2 + "' with linespoint title \"B\" lc rgb \"#009e73\" pt 2, '" + filename3 + "' with linespoint title \"T\" lc rgb \"#56b4e9\" pt 3";
		}
		else {
			message += ", '" + filename1 + "' with points title \"\" lc rgb \"dark-violet\" pt 1, '" + filename2 + "' with points title \"\" lc rgb \"#009e73\" pt 2, '" + filename3 + "' with points title \"\" lc rgb \"#56b4e9\" pt 3";
		}
	}
	message += ";";

	WriteToFile(message, "gnuplot.sh");

	return 0;
};