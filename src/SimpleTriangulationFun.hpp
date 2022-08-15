#ifndef SIMPLETRIANGULATION_FUN_HPP
#define SIMPLETRIANGULATION_FUN_HPP

#include <vector>
#include <string>
#include <fstream>

#include "SimpleTriangulationFun.hpp"


namespace SimpleTriangulation{

	void make_stl(Triangulation<std::vector<double>> t,std::string filename, std::string solidname);	
		
	Triangulation<std::vector<double>> PlanarTriangulation(
		std::vector<double> X,
		std::vector<double> Y,
		int N);

	void make_stl(Triangulation<std::vector<double>> t, std::string filename, std::string solidname = ""){
			/*
				Description: Constructs .stl ASCII format in filename from triangulation. If triangulation is 2D, an extra
					coordinate equal to zero is added to the vertices
				
				Input:
					- Triangulation<std::vector<double>> t
						t.faces.size() > 0
				
					- std::string filename
						data is saved in "filename", so .stl or other specifier must be included
					
					- std::string solidname
						name of solid at the top of the stl file
			
				Note: the triangulation is not necessarily oriented, so the unit normals might not belong to the
					same orientation. Make sure to call the orient() method before the make_stl() method if this behaviour
					is not desired
			*/
			
			auto exterior = [] (std::vector<double> a, std::vector<double> b){// exterior product lamda expression
			
				std::vector<double> n(3,0);
				n.at(0) = (a.at(1)*b.at(2) - a.at(2)*b.at(1)); 
				n.at(1) = (a.at(2)*b.at(0) - a.at(0)*b.at(2));
				n.at(2) = (a.at(0)*b.at(1) - a.at(1)*b.at(0));
				double n_mag = std::sqrt(n.at(0)*n.at(0) + n.at(1)*n.at(1) + n.at(2)*n.at(2));
				n.at(0) /= n_mag;
				n.at(1) /= n_mag;
				n.at(2) /= n_mag;
				
				return n;
			};
			
			auto subtract = [] (std::vector<double> a, std::vector<double> b){// lamda subtraction expression
				std::vector<double> result(a.size(),0);
				for (int i = 0; i < a.size(); i++){
					result.at(i) = a.at(i) - b.at(i);
				}
				return result;
			};
			
			/* If 2D converting to 3D */
			
			if (t.nodes.at(0).size() == 2){
				std::vector<std::vector<double>> new_nodes(t.nodes.size(),std::vector<double>(3,0));
				for (int i = 0; i< t.nodes.size(); i++){
					new_nodes.at(i).at(0) = t.nodes.at(i).at(0);
					new_nodes.at(i).at(1) = t.nodes.at(i).at(1);
					new_nodes.at(i).at(2) = 0;
				}
				t.nodes = new_nodes;
			}
			
			/* Writting to file */
			
			std::fstream file;
			file.open(filename,std::ios::out);
			
			std::vector<double> normal;
			std::vector<std::vector<double>> vertices;
			
			file << "solid " << solidname << std::endl;
			for (int i = 0; i < t.faces.size(); i++){
			
				vertices = t.get_vertices(i);
				
				normal = exterior(subtract(vertices.at(2),vertices.at(0)),subtract(vertices.at(1),vertices.at(0)));
				
				file << "facet normal " << normal.at(0) << " " << normal.at(1) << " " << normal.at(2) << std::endl;
				file << "outer loop" << std::endl;
					for (int j = 0; j < 3; j++){
						file << "vertex " << vertices.at(j).at(0) << " " << vertices.at(j).at(1) << " " << vertices.at(j).at(2) << std::endl;
					}
				file << "endloop" << std::endl;
				file << "endfacet" << std::endl;
			}
			file << "endsolid";
			file.close();	
	};
	
	Triangulation<std::vector<double>> PlanarTriangulation(
		std::vector<double> X,
		std::vector<double> Y,
		int N)
	{
		/*
			Description: Given the rectangle X \times Y \in R^2, it is subdivided into approximately N triangles
				which span the entire rectangle. First, each interval X,Y is subdivided into M = ceiling(sqrt(N/2)) 
				subintervals. These subintervals partition the rectangle into M^2 subrectangles. Each subrectangle 
				is divided into two triangles by joining opposite vertices. 
				
				Conclusively, given N, a triangulation of X \times Y into 2*ceil(sqrt(N/2))^2 triangles is produced
			
			Input:
			
				- std::vector<double> X
					First side of rectangle: [X.at(0), X.at(1)]
				- std::vector<double> Y
					Second side of rectangle: [Y.at(0), Y.at(1)]
				- int N
					A triangulation consisting of 2*ceil(sqrt(N/2))^2 rectangles is produced
			
			Output:
			
				- std::vector<std::vector<std::vector<double>>> triangles
					triangles.at(i) is the ith triangle consisting of vertices triangles.at(i).at(0), triangles.at(i).at(1),
					triangles.at(i).at(2)
			
		*/
		
		/* Subdividing X,Y */
		
		int M = std::ceil(std::sqrt(N/2)); //number of subintervals in X,Y
		
		std::vector<double> X_discrete(M+1,0);
		std::vector<double> Y_discrete(M+1,0);
		
		for (int i = 0; i < (M+1); i++){
			X_discrete.at(i) = X.at(0) + (X.at(1)-X.at(0))*i/M;
			Y_discrete.at(i) = Y.at(0) + (Y.at(1)-Y.at(0))*i/M;
		}
		
		/* Constructing Triangles */
		
		std::vector<std::vector<std::vector<double>>> triangles(2*M*M,std::vector<std::vector<double>>(3,std::vector<double>(2,0))); //triangles initialized to all zeros
		
		int count = 0;// counts rectangles
		for (int i = 0; i < M; i++){
			for (int j = 0; j < M; j++){
				// (i,j) denotes the (i,j)th discrete subrectangle of X \times Y
				
				/* lower right triangle */
				
				triangles.at(count).at(0) = std::vector<double> {X_discrete.at(i), Y_discrete.at(j)};
				triangles.at(count).at(1) = std::vector<double> {X_discrete.at(i+1), Y_discrete.at(j)};
				triangles.at(count).at(2) = std::vector<double> {X_discrete.at(i+1), Y_discrete.at(j+1)};
				
				count++;
				
				/* upper left triangle */
				
				triangles.at(count).at(0) = std::vector<double> {X_discrete.at(i), Y_discrete.at(j)};
				triangles.at(count).at(1) = std::vector<double> {X_discrete.at(i), Y_discrete.at(j+1)};
				triangles.at(count).at(2) = std::vector<double> {X_discrete.at(i+1), Y_discrete.at(j+1)};
				
				count++;
			}
		}
		
		return(Triangulation<std::vector<double>>());
		
	}// PlanarTriangulation()
	
	
}// namespace

#endif // SIMPLETRIANGULATION_FUN_HPP
