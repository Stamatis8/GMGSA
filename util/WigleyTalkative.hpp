#include <string>
#include <string>
#include "../modelers/WigleyModeler.hpp"

std::vector<std::string> WigleyTalkative(WigleyModeler m, std::string options) {
/*
	Description: Constructs titles/ parameter labels/ etc for graphs involving the WigleyModeler modeler

	Input:
		- std::string options
			options == "title"
				Construct a title whose first line mentions the parameters, then ranges for each parameter
				with at most three ranges per line and finally constant values of parameters with at most three per line.
				This string is contained in output.at(0)
			options == "labels"
				output.at(i) = ith parameter label
		
*/

	std::vector<std::string> output;
	auto D = m.design_space();// design space
	auto P = m.get_particulars();// {L,B,d,c1,c2,c3}

	if (options == "title") {
		output = { "" };

		if (m.type == 0) {
			output.at(0) = "Wigley Hull with L,B,T parameters\\n";
			output.at(0) += " L in [" + std::to_string(D.at(0).at(0)) + "," + std::to_string(D.at(0).at(1)) + "], B in [" + std::to_string(D.at(1).at(0)) + "," + std::to_string(D.at(1).at(1)) + "], T in [" + std::to_string(D.at(2).at(0)) + "," + std::to_string(D.at(2).at(1)) + "]\\n";
			output.at(0) += " c1 = " + std::to_string(P.at(3)) + ", c2 = " + std::to_string(P.at(4)) + ", c3 = " + std::to_string(P.at(5)) + "\\n";
		}
		else if (m.type == 1) {
			output.at(0) = "Wigley Hull with L,B,T,c1,c2,c3 parameters\\n";
			output.at(0) += " L in [" + std::to_string(D.at(0).at(0)) + "," + std::to_string(D.at(0).at(1)) + "], B in [" + std::to_string(D.at(1).at(0)) + "," + std::to_string(D.at(1).at(1)) + "], T in [" + std::to_string(D.at(2).at(0)) + "," + std::to_string(D.at(2).at(1)) + "]\\n";
			output.at(0) += " c1 in [" + std::to_string(D.at(3).at(0)) + "," + std::to_string(D.at(3).at(1)) + "], c2 in [" + std::to_string(D.at(4).at(0)) + "," + std::to_string(D.at(4).at(1)) + "], c3 in [" + std::to_string(D.at(5).at(0)) + "," + std::to_string(D.at(5).at(1)) + "]\\n";
		}
		else if (m.type == 2) {
			output.at(0) = "Wigley Hull with B/L, T/L parameters\\n";
			output.at(0) += " B/L in [" + std::to_string(D.at(0).at(0)) + "," + std::to_string(D.at(0).at(1)) + "], T/L in [" + std::to_string(D.at(1).at(0)) + "," + std::to_string(D.at(1).at(1)) + "]\\n";
			output.at(0) += " L = " + std::to_string(P.at(0)) + ", c1 = " + std::to_string(P.at(3)) + ", c2 = " + std::to_string(P.at(4)) + ", c3 = " + std::to_string(P.at(5)) + "\\n";
		}
		else if (m.type == 3) {
			output.at(0) = "Wigley Hull with B/L, T/L, c1, c2, c3 parameters\\n";
			output.at(0) += " B/L in [" + std::to_string(D.at(0).at(0)) + "," + std::to_string(D.at(0).at(1)) + "], T/L in [" + std::to_string(D.at(1).at(0)) + "," + std::to_string(D.at(1).at(1)) + "]\\n";
			output.at(0) += " c1 in [" + std::to_string(D.at(3).at(0)) + "," + std::to_string(D.at(3).at(1)) + "], c2 in [" + std::to_string(D.at(4).at(0)) + "," + std::to_string(D.at(4).at(1)) + "], c3 in [" + std::to_string(D.at(5).at(0)) + "," + std::to_string(D.at(5).at(1)) + "]\\n";
			output.at(0) += " L = " + std::to_string(P.at(0)) + "\\n";
		}
	}
	else if (options == "labels") {
		if (m.type == 0) {
			output.push_back("L");output.push_back("B");output.push_back("T");
		}
		else if (m.type == 1) {
			output.push_back("L");output.push_back("B");output.push_back("T");output.push_back("c_1");output.push_back("c_2");output.push_back("c_3");
		}
		else if (m.type == 2) {
			output.push_back("B/L");output.push_back("T/L");
		}
		else if (m.type == 3) {
			output.push_back("B/L");output.push_back("T/L");output.push_back("c_1");output.push_back("c_2");output.push_back("c_3");
		}
	}

	return output;
}