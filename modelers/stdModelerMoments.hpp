#ifndef STDMODELERMOMENTS_HPP
#define STDMODELERMOMENTS_HPP

#include <vector>

#include "../src/geom_moments/geom_moments.hpp"
#include "../src/smpl_triangulation/smpl_triangulation.hpp"

template<typename PM>
class stdModelerMoments: public PM{
	/*
		Description: This class is an extension for some arbitrary base parametric modeler class PM. This class allows for
			the calculation of geometric moments of PM, via first triangulating PM and then approximating the requested moment
			for the current design in PM.
			
		PM typename requirements:
		
			PM.evaluate(std::vector<double> args)
				- returns an std::vector<double> with evaluated points for current design
			
			PM.domain()
				- retrurns std::vector<std::vector<double>> with each element being the upper and lower bound for each
					surface/design parameter for the current design.
				- domain must have two elements (ie a surface design)
			
		Notes:
			
			- Moments are evaluated with geom_moments library
			
			- Triangulation is constructed with smpl_triangulation library
			
	*/
public:

	/* Constructors */

	stdModelerMoments(PM modeler, int N_triangles_in = 100): PM(modeler){
		/*
			Description: Constructor which passes modeler argument to default copy-constructor of PM
		*/
		
		this->N_triangles = N_triangles_in;
	}
	
	double moment(int p, int q, int r, bool is_translation_invariant = false, bool is_scaling_invariant = false){
		/*
			Description:
				  - Calculates the s = p + q + r order geometric moment of the current design in modeler.
				  - is_translation_invariant  == true calculates the translation invariant of said moment
				  - is_scaling_invariant == true calculates scaling invariant of said moment
		*/
		
		return MomentSthOrder(this->triangulate(),p,q,r,p+q+r,is_translation_invariant,is_scaling_invariant);
	}
	
	std::vector<std::vector<std::vector<double>>> triangulate(std::string filename = "", std::string solidname = ""){
		/*
			Description: Construct a triangulation for this->design in approximately N triangles. If filename != "" constructs .stl file
				filename
				
			Input:
				- filename
					if filename != "", makes filename in STL ASCII format. Else, no STL file is created
					".stl" ending must be included, ie filename = "myfile.stl"
					
				- solidname
					if STL file is to be constructed, this will be the solid name
					
			Output:
				- std::vector<std::vector<std::vector<double>>> T
					T.at(i) is the ith triangle with its 3 elements being its three vertices
						
		*/
		
		/* Triangulate domain */
		
		std::vector<std::vector<double>> d = this->domain();
		
		smpl_triangulation::Triangulation<std::vector<double>> T = smpl_triangulation::PlanarTriangulation(d.at(0),d.at(1), this->N_triangles);
	
		/* mapping 2D triangulation on design */
		
		for (int vertex = 0; vertex < T.nodes.size(); vertex++){
			T.nodes.at(vertex) = this->evaluate({T.nodes.at(vertex).at(0),T.nodes.at(vertex).at(1)});//domain to design
		}
		
		if (filename != ""){// make STL file
			smpl_triangulation::make_stl(T,filename,solidname);
		}
		
		return T;// cast of T to vec<vec<vec>>> format
	}
	
private:

	int N_triangles;// Approx. number of triangles to construct modeler mesh
	
};
#endif //STDMODELERMOMENTS_HPP
