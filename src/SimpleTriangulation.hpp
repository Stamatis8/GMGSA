#ifndef SIMPLETRIANGULATION_HPP
#define SIMPLETRIANGULATION_HPP

#include <vector>
#include <cmath>

std::vector<std::vector<std::vector<double>>> SimpleTriangulation(
	std::vector<double> X,
	std::vector<double> Y,
	int N);

std::vector<std::vector<std::vector<double>>> SimpleTriangulation(
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
	
	return(triangles);
	
}// SimpleTriangulation()

#endif// SIMPLETRIANGULATION_HPP
