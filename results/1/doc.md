27/8/22
# Description
The ith dat file contains two columns. The first is the number of samples and the second is the sensitivity index for the ith parameter.

	WigleyModeler Wigley {{0.8,1.2},{0.08,0.12},{0.05,0.075},0.2,0,1};

	WigleyAnalyticMoments WigleyAnalytic {Wigley};

GMGSA():

	int sub_population_size = 2;// see references in DPS.hpp
	int max_iterations = 3;// see references in DPS.hpp
	double omega = 1;// see references in DPS.hpp

DPS() was used for sampling generation

order of SSVs was 4