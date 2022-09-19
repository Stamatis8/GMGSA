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
//#define SSV_NONE_SCALING_INVARIANT

#include "../util/timer.hpp"
#include "../util/WriteToFile.hpp"
#include "../util/WigleyTalkative.hpp"

#include "../src/GMGSA.hpp"
#include "../modelers/WigleyModeler.hpp"
#include "../modelers/WigleyAnalyticMoments.hpp"

int main() {

	/* Experiment initialization */

	std::vector<long int> N_range = { 2000, 20000 };// SSV order range of experiment
	long int N_step = 2000;// calculate at N_range.at(0) + i*N_step until N_range.at(1)
	int order = 2;// number of samples to be used in experiment
	int runs = 1;// number of runs of experiment

	std::srand(11);// random seed
	//2321

	/* Wigley modeler initialization */

	//old
	//WigleyModeler model{ { 0.8,1.2 }, { 0.08,0.12 }, { 0.033,0.046 }, 0.2, 0,  1 };

	//WigleyModeler model{ { 0.8,1.2 }, { 0.08,0.12 }, { 0.033,0.046 }, 0.2, 0,  1 };

	double c = 1;
	WigleyModeler model{ { c*0.8,c*1.2 }, { c*0.08,c*0.12 }, { c*0.033,c*0.046 }, 0.2, 0,  1 };

	//WigleyModeler model{ 0.1,{0.06,0.15},{0.06,0.15},0.2,0,1,"ratios" };

	//L,B,T
	//WigleyModeler model{ { 0.8,1.2 }, { 0.072,0.12 }, { 0.06,0.216 }, 0.2, 0,  1 };

	int dim = model.design_space().size();//number of dimensions

	/* Graph title and label creation */

	std::string title, xlabel, ylabel;

	title = WigleyTalkative(model, "title").at(0);

	title += "using simple MonteCarlo sampling. SSV order = " + std::to_string(order) + "\\n";
	xlabel = "order";
	ylabel = "SI";


	/* Wigley moment calculator initialization */

#ifdef GMGSA_USE_BOOST_MULTIPRECISION
	typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float< 30 > > BigFloat;
#else
	typedef long double BigFloat;
#endif

	WigleyAnalyticMoments<BigFloat> Wigley{ model };

	/* Calculation */

	std::vector<std::vector<std::vector<BigFloat>>> P_data(dim);// P_data.at(i) = ith parameter data
	std::vector<BigFloat> SI;//Sensitivity indices

	timer t;
	for (int run = 1; run <= runs; run++) {
		std::cout << "Starting run (" << std::to_string(run) << " / " << std::to_string(runs) << ")" << std::endl;

		for (int N = N_range.at(0); N <= N_range.at(1); N += N_step) {
			t.begin();

			SI = GMGSA<WigleyAnalyticMoments<BigFloat>, BigFloat>(Wigley, N, order);
			
			for (int par = 0; par < P_data.size(); par++) {
				P_data.at(par).push_back({ double(N),SI.at(par),double(order) });
			}

			std::cout << "Finished with N = " << N << " in ";
			t.display();
			std::cout << "s " << std::endl;
		}

		for (int par = 0; par < P_data.size(); par++) {
			WriteToFile(P_data.at(par), "data/Par" + std::to_string(par + 1) + "_" + std::to_string(run) + ".dat");
			P_data.at(par).clear();
		}

		std::cout << "Finished run (" << std::to_string(run) << " / " << std::to_string(runs) << ")" << std::endl << std::endl;
	}


	/* Constructing gnuplot script */

	std::string message =
		"set title \"" + title + "\"; set xlabel \"" + xlabel + "\"; set ylabel \"" + ylabel + "\" rotate by 0;"
		"set key outside;"
		;

	std::string filename;
	std::vector<std::string> labels = WigleyTalkative(model, "labels");
	std::vector<std::string> colours = { "dark-violet","#009e73","#56b4e9","#e69f00","#f0e442","#0072b2" };
	for (int run = 1; run <= runs; run++) {

		if (run == 1) {
			for (int par = 0; par < dim; par++) {

				filename = "data/Par" + std::to_string(par + 1) + "_" + std::to_string(run) + ".dat";

				if (par == 0) {
					message += "plot '" + filename + "' with linespoint title \"" + labels.at(par) + "\" lc rgb \"" + colours.at(par) + "\" pt " + std::to_string(par + 1);
				}
				else {
					message += ", '" + filename + "' with linespoint title \"" + labels.at(par) + "\" lc rgb \"" + colours.at(par) + "\" pt " + std::to_string(par + 1);
				}
			}
		}
		else {
			for (int par = 0; par < dim; par++) {

				filename = "data/Par" + std::to_string(par + 1) + "_" + std::to_string(run) + ".dat";

				message += ", '" + filename + "' with points title \"\" lc rgb \"" + colours.at(par) + "\" pt " + std::to_string(par + 1);

			}
		}
	}
	message += ";";

	WriteToFile(message, "gnuplot.sh");

	return 0;
};