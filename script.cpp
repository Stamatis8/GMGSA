#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include "src/smpl_triangulation/smpl_triangulation.hpp"
#include "src/geom_moments/geom_moments.hpp"
#include "src/geom_moments/stl_util/stl2vec.hpp"
#include "util/WriteToFile.hpp"

#include "src/DPS.hpp"
#include "src/Wigley_Modeler.hpp"
#include "src/GMGSA.hpp"

using namespace smpl_triangulation;

int main() {

	int test = 12;
	if (test == 1){
		// DPS test
			
		std::vector<std::vector<double>> X = {{0,1},{0,1},{0,1}};
			
		std::vector<std::vector<double>> S = DPS(X,std::vector<std::vector<double>>(),20,1,5,1);
			
		WriteToFile(S,"samples.dat");	
	}
	else if (test == 12){
	WigleyModeler model(100,40,10);
	int N = 100;
	int N_t = 50;
	int order = 0;
	std::vector<double> SI = GMGSA<WigleyModeler>(model,N,N_t,order);
	
	std::cout << "Order " << order << " SI: " << std::endl;
		for (int i = 0; i < SI.size(); i++){
			std::cout << SI.at(i) << std::endl;
		}
	}
	else if (test == 11){
		WigleyModeler model(100,30,10);
		int N_triangle = 1000;
		std::vector<std::vector<std::vector<double>>> triangles = model.triangulate(N_triangle,"Wigley" + std::to_string(N_triangle) +".stl","Wigley");
	
		double V = MomentSthOrder(triangles,0,0,0,0,false,false);
		double Cx = MomentSthOrder(triangles,1,0,0,1,false,false)/V;
		double Cy =	MomentSthOrder(triangles,0,1,0,1,false,false)/V;
		double Cz = MomentSthOrder(triangles,0,0,1,1,false,false)/V;
		
		std::cout << "Volume is: "<< V <<std::endl;
		std::cout << "Centoid is at: ( "<<Cx<<" , "<<Cy<<" , "<<Cz<<" )"<<std::endl;
	}
	else if (test == 10){
		std::vector<double> A {1,2,3,4,5,6,7};
		std::vector<double> a = A;
		a.erase(a.begin()+2);
		
		std::cout << "A: \n";
		for (int i = 0; i< A.size(); i++){
			std::cout << A.at(i) << std::endl;
		}
		
		std::cout << "a: \n";
		for (int i = 0; i< a.size(); i++){
			std::cout << a.at(i) << std::endl;
		}
	}
	else if (test == 9){
		int order = 4;
		std::string filename = "botton1.stl";
		std::vector<double> V = SSV(filename,order);
		
		std::cout << "Order " << order << " SSV of " << filename << "is: " << std::endl;
		for (int i = 0; i < V.size(); i++){
			std::cout << V.at(i) << std::endl;
		}
	}
	else if (test == 8){
	
		std::string filename = "botton1.stl";
	
		double V = MomentSthOrder(filename,0,0,0,0,false,false);
		double Cx = MomentSthOrder(filename,1,0,0,1,false,false)/V;
		double Cy =	MomentSthOrder(filename,0,1,0,1,false,false)/V;
		double Cz = MomentSthOrder(filename,0,0,1,1,false,false)/V;
		
		std::cout << "Volume of .stl file is: "<< V <<std::endl;
		std::cout << "Centoid is at: ( "<<Cx<<" , "<<Cy<<" , "<<Cz<<" )"<<std::endl;
	}
	else if (test == 2){
		// SimpleTriangulation test
	
		/*
		std::vector<std::vector<std::vector<double>>> triangles = SimpleTriangulation(std::vector<double> {0, 1}, std::vector<double> {0, 1}, 100);
		
		std::cout << std::endl;
		for (int i = 0; i < triangles.size(); i++){
			std::cout << "Triangle number " + std::to_string(i+1) + " :\n";
			for (int j = 0; j < 3; j++){
			std::cout << "( " + std::to_string(triangles.at(i).at(j).at(0)) + " ,"
						  	+ std::to_string(triangles.at(i).at(j).at(1))
						  	+ " )" << std::endl;
			}
			std::cout << std::endl << std::endl;
		}
		std::cout<< std::endl;
		*/
	}
	else if (test == 3){
		// NchooseK struct test
		
		NchooseK_cache NchooseK;
		
		std::cout << "4 choose 2 is : " << NchooseK.get(4,2) << std::endl;
		
		WriteToFile(NchooseK.nk,"NchooseK_4_2.txt");
		
		NchooseK.get(16,1);
		
		WriteToFile(NchooseK.nk,"NchooseK_16_1.txt");
		
		NchooseK.get(3,1);
		
		std::cout<<"Size of nk after (16,1) (1,3) commands: "<<NchooseK.nk.size()<<std::endl;
		
		NchooseK.get(1000,1);
		
		WriteToFile(NchooseK.nk,"NchooseK_large.txt");
	
		std::cout<<"1521 choose 3 is equal to: "<<NchooseK.get(1521,3)<<std::endl;
	}
	else if (test == 4){
		// MomentSthOrder test
		
		std::vector<std::vector<std::vector<double>>> triangles = {
																	{//triangle 1
																	  {0,0,0
																	  },
																	  {1,0,0
																	  },
																	  {1,1,0
																	  }
																	},
																	{//triangle 2
																	  {0,0,0
																	  },
																	  {1,1,0
																	  },
																	  {0,1,0
																	  }
																	},
																	{//triangle 3
																	  {0,0,0
																	  },
																	  {1,0,1
																	  },
																	  {1,0,0
																	  }
																	},
																	{//triangle 4
																	  {0,0,0
																	  },
																	  {0,0,1
																	  },
																	  {1,0,1
																	  }
																	},
																	{//triangle 5
																	  {1,0,0
																	  },
																	  {1,0,1
																	  },
																	  {1,1,0
																	  }
																	},
																	{//triangle 6
																	  {0,0,0
																	  },
																	  {0,1,0
																	  },
																	  {0,0,1
																	  }
																	},
																	{//triangle 7
																	  {1,0,1
																	  },
																	  {0,0,1
																	  },
																	  {1,1,0
																	  }
																	},
																	{//triangle 8
																	  {0,0,1
																	  },
																	  {0,1,0
																	  },
																	  {1,1,0
																	  }
																	}
																  };
		/*
		triangles = {
						{//triangle 1
							{0,0,0
							},
							{1,0,0
							},
							{0,1,0
							}
						},
						{//triangle 2
							{0,0,0
							},
							{0,1,0
							},
							{0,0,1
							}
						},
						{//triangle 3
							{0,0,0
							},
							{0,0,1
							},
							{1,0,0
							}
						},
						{//triangle 4
							{0,0,1
							},
							{0,1,0
							},
							{1,0,0
							}
						}
					};											  
		*/
		
		std::cout<< "With New Moment Function: " << MomentSthOrder(triangles,0,0,0,0,false,false)<<std::endl;
	}
	else if (test == 5){
		
		std::vector<std::vector<double>> nodes = {
													{0,0,0},
													{1,0,0},
													{0,1,0},
													{1,1,0},
													{2,0,0},
													{2,1,0}
											};
		std::vector<std::vector<int>> faces = {
												{0,2,1,1,0,0},
												{1,3,2,1,0,2},
												{1,3,5,2,3,1},
												{4,5,1,2,3,3}
											};
		Triangulation<std::vector<double>> T(nodes,faces);
		
		T.orient();
		
		make_stl(T,"mytriangulation.stl");
		
		int temp = 10;
	}
	else if (test == 6){
		Triangulation<std::vector<double>> T = PlanarTriangulation({0,2},{0,2},8);
		
		T.orient();
		
		make_stl(T,"mytriangulation.stl");
	}
	else if (test == 7){
		J_cache J;
		
		std::cout << "J(0,0) = "<< J.get(0,0) << std::endl;
		std::cout << "J(1,0) = "<< J.get(1,0) << std::endl;
		std::cout << "J(0,1) = "<< J.get(0,1) << std::endl;
		std::cout << "J(1,1) = "<< J.get(1,1) << std::endl;
		std::cout << "J(2,0) = "<< J.get(2,0) << std::endl;
		std::cout << "J(0,2) = "<< J.get(0,2) << std::endl;
	}
	return 0;
}
