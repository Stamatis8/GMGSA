#include <iostream>
#include <vector>

#include "src/DPS.hpp"
#include "src/SimpleTriangulation.hpp"
#include "src/MomentSthOrder.hpp"

#include "util/WriteToFile.hpp"
#include "util/NchooseK.hpp"


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
																	  {0,1,0
																	  },
																	  {1,1,0
																	  }
																	},
																	{//triangle 3
																	  {0,0,0
																	  },
																	  {1,0,0
																	  },
																	  {1,0,1
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
																	  {0,0,0
																	  },
																	  {1,0,0
																	  },
																	  {0,0,1
																	  }
																	},
																	{//triangle 6
																	  {0,1,0
																	  },
																	  {1,1,0
																	  },
																	  {0,1,1
																	  }
																	},
																	{//triangle 7
																	  {0,0,1
																	  },
																	  {1,0,0
																	  },
																	  {1,1,0
																	  }
																	},
																	{//triangle 8
																	  {0,0,1
																	  },
																	  {0,1,1
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
			
		std::cout<< MomentSthOrder(triangles,0,0,0,false)<<std::endl;
	}
	return 0;
}
