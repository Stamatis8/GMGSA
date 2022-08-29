#ifndef WRITETOFILE_HPP
#define WRITETOFILE_HPP

#include <string>
#include <vector>
#include <fstream>

template<typename scalar>
bool WriteToFile(
	std::vector<std::vector<scalar>> X,
	std::string filename,
	std::string message = ""
	);

template<typename scalar>
bool WriteToFile(
	std::vector<std::vector<scalar>> X,
	std::string filename,
	std::string message
	)
{
	/*
		Description: Writes each element of X to a row in filename. Each element of X.at(i) is separated by
			a space in the respective line in filename
			
		Input:
			
			- std::vector<std::vector<double>> X
				X is written in filename. X.at(i) is written at ith row of filename. Each element of X.at(i) is
					separated by a space
			- std::string filename
				file to write X at. Must include extension (ie "data.txt" or "data.dat")
	
		Output:
		
			- bool out
				Signifies completion
				
		Notes: If filename exists, then it is overwritten
	*/
	
	std::fstream file;
	bool out;
	
	file.open(filename,std::ios::out);
	if (file.is_open()){
	
		out = true;
	
		if (message != "") {//print message
			file << "#" << message << std::endl;
		}

		for (int i = 0;i < X.size();i++){
			for (int j = 0; j < X.at(i).size(); j++){
				file << X.at(i).at(j);
				if (j != (X.at(i).size() - 1)) file << " ";// Do not add space after last element of each row in file
			}
			if (i != (X.size() - 1)) file << std::endl;// Do not add newline after last row in file
		}
	}
	else{
		out = false;	
	}
	file.close();	
	
	return(out);
	
}// WriteToFile()

#endif// WRITETOFILE_HPP
