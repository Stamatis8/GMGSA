#include <iostream>
#include <vector>
#include <cmath>

#include "src/DPS.hpp"
#include "src/SimpleTriangulationClass.hpp"
#include "src/SimpleTriangulationFun.hpp"
#include "src/MomentSthOrder.hpp"

#include "util/WriteToFile.hpp"
#include "util/NchooseK.hpp"

using namespace SimpleTriangulation;

int main() {
	
	int test = 4;
	if (test == 1){
		// DPS test
			
		std::vector<std::vector<double>> X = {{0,1},{0,1},{0,1}};
			
		std::vector<std::vector<double>> S = DPS(X,std::vector<std::vector<double>>(),20,1,5,1);
			
		WriteToFile(S,"samples.dat");	
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
		
		
			
		std::cout<< "With Moment Function: " << old_MomentSthOrder(triangles,0,0,0,false)<<std::endl;
		
		/* DEBUG: Manual Calculation */
		
		double M = 0;
		for (int t = 0; t < triangles.size(); t++){
			double I_t = 0;
			
			std::vector<double> a(3,0);
			std::vector<double> b(3,0);
			std::vector<double> c(3,0);
			std::vector<double> n(3,0);
			
			double a_mag;
			double b_mag;
			double n_mag;
			
			a.at(0) = triangles.at(t).at(1).at(0) - triangles.at(t).at(0).at(0);
			a.at(1) = triangles.at(t).at(1).at(1) - triangles.at(t).at(0).at(1);
			a.at(2) = triangles.at(t).at(1).at(2) - triangles.at(t).at(0).at(2);
			
			b.at(0) = triangles.at(t).at(2).at(0) - triangles.at(t).at(0).at(0);
			b.at(1) = triangles.at(t).at(2).at(1) - triangles.at(t).at(0).at(1);
			b.at(2) = triangles.at(t).at(2).at(2) - triangles.at(t).at(0).at(2);
			
			c = triangles.at(t).at(0); 
			
			a_mag = std::sqrt(a.at(0)*a.at(0) + a.at(1)*a.at(1) + a.at(2)*a.at(2));
			b_mag = std::sqrt(b.at(0)*b.at(0) + b.at(1)*b.at(1) + b.at(2)*b.at(2));
				
			n.at(0) = (a.at(1)*b.at(2) - a.at(2)*b.at(1)); 
			n.at(1) = (a.at(2)*b.at(0) - a.at(0)*b.at(2));
			n.at(2) = (a.at(0)*b.at(1) - a.at(1)*b.at(0));
			n_mag = std::sqrt(n.at(0)*n.at(0) + n.at(1)*n.at(1) + n.at(2)*n.at(2));
			n.at(0) /= n_mag;
			n.at(1) /= n_mag;
			n.at(2) /= n_mag;
			
			double ab = a.at(0)*b.at(0) + a.at(1)*b.at(1) + a.at(2)*b.at(2);
			double VDa = std::sqrt(a_mag*a_mag*b_mag*b_mag - ab*ab);
			
			I_t = (1/6)*(n.at(0)*(a.at(0)/3+b.at(0)/3+c.at(0)) 
					   + n.at(1)*(a.at(1)/3+b.at(1)/3+c.at(1))
					   + n.at(2)*(a.at(2)/3+b.at(2)/3+c.at(2))	
						);
			
			I_t *= VDa;
			
			M += I_t;
		}	
		
		std::cout<< "With Manual Calculation: " << M <<std::endl;
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
	return 0;
}
