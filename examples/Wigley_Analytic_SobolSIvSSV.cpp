/*
	Description: Perform Sensitivity analysis for the Wigley hull for a fixed number of samples and SSV in a certain range.
		Capability for multiple runs.
*/


#include <vector>
#include <iostream>
#include <iomanip>
#include <string>

#include "../util/timer.hpp"
#include "../util/WriteToFile.hpp"
#include "../util/WigleyTalkative.hpp"

#include "../modelers/WigleyModeler.hpp"
#include "../modelers/WigleyAnalyticSobol.hpp"



int main() {

	/* Experiment initialization */

	double order_max = 15;// SI will be calculated for SSV orders 0,...,order_max

	std::srand(11);// random seed 
	//68
	//2321

	/* Wigley modeler initialization */

	//L,B,T
	WigleyModeler model{ { 0.75,1.25 }, { 0.072,0.12 }, { 0.103,0.1725 }, 0.2, 0,  1 };

	int dim = model.design_space().size();//number of dimensions

	/* Graph title and label creation */

	std::string title, xlabel, ylabel;

	title = WigleyTalkative(model, "title").at(0);
	
	title += "using analytic Sobol indices\\n";
	xlabel = "order";
	ylabel = "SI";

	/* Calculation */
    
    typedef long double BigFloat;
	
    std::vector<std::vector<std::vector<BigFloat>>>P_data(dim);// P_data.at(i) = ith parameter data
	std::vector<BigFloat> SI;//Sensitivity indices
    std::vector<std::vector<BigFloat>> Analytic_SI = WigleyAnalyticSobol<BigFloat>(model,order_max);

    for (int order = 0; order <= order_max; order++) {
        SI = Analytic_SI.at(order);

        for (int par = 0; par < P_data.size(); par++) {
            P_data.at(par).push_back({ double(order),SI.at(par)});
        }
    }

    for (int par = 0; par < P_data.size(); par++) {
        WriteToFile(P_data.at(par), "data/Par" + std::to_string(par + 1) + ".dat");
        P_data.at(par).clear();
    }


	/* Constructing gnuplot script */

	std::string message =
		"set title \"" + title + "\"; set xlabel \"" + xlabel + "\"; set ylabel \"" + ylabel + "\" rotate by 0;"
		"set key outside;"
		;

	std::string filename;
	std::vector<std::string> labels = WigleyTalkative(model, "labels");
	std::vector<std::string> colours = { "dark-violet","#009e73","#56b4e9","#e69f00","#f0e442","#0072b2" };


    for (int par = 0; par < dim; par++) {

        filename = "data/Par" + std::to_string(par+1) + ".dat";

        if (par == 0) {
            message += "plot '" + filename + "' with linespoint title \"" + labels.at(par) + "\" lc rgb \"" + colours.at(par) + "\" pt " + std::to_string(par + 1);
        }
        else {
            message += ", '" + filename + "' with linespoint title \"" + labels.at(par) + "\" lc rgb \"" + colours.at(par) + "\" pt " + std::to_string(par + 1);
        }
    }

	
	message += ";";

	WriteToFile(message, "gnuplot.sh");

	return 0;
};