#include <tinynurbs/tinynurbs.h>
#include <glm/glm.hpp>
#include <iostream>

int main() {
	
	tinynurbs::Curve<float> crv; // Planar curve using float32
	crv.control_points = {glm::vec3(-1, 0, 0), // std::vector of 3D points
                      	glm::vec3(0, 1, 0),
                      	glm::vec3(1, 0, 0)
                     	};
	crv.knots = {0, 0, 0, 1, 1, 1}; // std::vector of floats
	crv.degree = 2;

	std::cout << "Hello World!" << std::endl;
	
	return 1;
}
