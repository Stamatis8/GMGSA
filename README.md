# Geometric Moment-dependent Global Sensitivity Analysis (GMGSA)

- C++ header library to preform geometric-moment based GSA as in [1]
- This library attaches to an existant parametric modeler, for which the requirements are listed below
- The parametric modeler is used as `class template` with minimal requirements, which means that the wrapper around already existing modelers should be is easy to build

# Dependancies

- Standard Library
- This [https://github.com/Stamatis8/geom_moments](https://github.com/Stamatis8/geom_moments) repository. Used to calculate geometric moments. It is already included in `src` file

### Optional

- This [https://github.com/Stamatis8/smpl_triangulation](https://github.com/Stamatis8/smpl_triangulation) repository. Used in examples to constuct triangulation of surface from parametrized model. It is also used to construct .STL files from the various designs. It is already included in `src` file

# Usage

All files that are needed to use this library are located in `src`. Move them at some `path` where the compiler has access and then simply `#include "path/GMGSA.hpp"` at the file you want to use this package from (ie copy-paste all files from `src` in your project folder and then `#include "GMGSA.hpp"`).

# Overview

This library is built around the `GMGSA()` function which given a parametric modeller returns the generalized total sensitivity index of each parameter. To describe the various components, we will follow a top-down approach. Having read the *Description of GMGSA approach* section (see below), would be useful. The generalized total sensitivity indices for each parameter are approximated using [4], in `GSI_T_estimator.hpp`. The inputs to this algorithm are the Shape Signature Vectors (SSV) for a number of randomly sampled designs from the design space. To calculate SSV, the `geom_moments` library is used, which given a triangulated mesh provides SSV of arbitrary order. For the sampling, the _Dynamic Propagation Sampling_ technique [2] is used which is implemented in `DPS.hpp`. Finally, the task of translating a given model to a triangulated mesh is left up to the user. However, a suggested solution is provided and demonstrated in the Wigley Hull example, which uses the `smpl_triangulation` library. 

We will now discuss the Wigley hull example

## Wigley hull example

The modifed Wigley hull-form is well-known and widely used for experimental and numerical studies [6]. The goal is to perform GMGSA of order 2 on a Wigley hull parametric modeler. The proposed parametric modeler for the wigley hull can be found in `modelers/WigleyModeler.hpp`. The example can be found in `examples/Wigley.cpp`.

	#include <vector>
	#include <iostream>
	
	#include "../src/GMGSA.hpp"
	
	#include "../modelers/WigleyModeler.hpp"
	
	int main(){
	
		double L = 100;// length
		double B = 30;// breadth
		double d = 15;// depth
	
		WigleyModeler model(L,B,d);
		
		std::vector<double> GSI; //generalized sensitivity indices
		
		GSI = GMGSA(model,50,200,2);// 50 samples, 200 triangles per sample, order 2 SSV
		
		for (int i = 0; i < GSI.size(); i++){	
			std::cout << "Sensitivity index for t" << i << " is: " << GSI.at(i) << std::endl;
		}
	
		return 0;
	}

# Documentation

Unfortunately, no organized documentation document exists yet. However each object is documented thoroughly at its definition. What follows is a listing of the requirements for the parametric modeler class

## Parametric modeler requirements

The `GMGSA()` function has the following prototype:

	template<class PM>
	std::vector<double> GMGSA(PM modeler, int N ,int N_triangles, int order);
	
`N` is the number of random samples to generate, `N_triangles` the number of triangles in each mesh and `order`, the order of the SSV to use. Now, the `modeler` variable is of type `PM`. The requirements for `PM` are listed in `src/GMGSA.hpp` as follows:
	
	class modeler requirements:
		
			- modeler.design_space()
				- returns an std::vector<std::vector<double>> vec
				- vec.at(i) is of size 2 and is equal to the range of the ith parameter
				- ie t_i \in [vec.at(i).at(0), vec.at(i).at(1)]
				- therefore it must also be true that number of parameters == vec.size()
				- modeler.design_space().size() > 0
			
			- modeler.set_design(std::vector<double> design)
				- design.size() == modeler.design_space().size()
				- design contains a set of parameters for the modeler
				- the modeler class must save this design until it is changed again
				
			- modeler.triangulate(N_triangles)
				- returns a triangulation of the currently set design in *approximately* N_triangles
				- the triangulation is returned in std::vector<std::vector<std::vector<double>>> triang format
				- the triangulation must be oriented
				- triang.at(i) is the ith triangle
				- triang.at(i).at(j) is the jth vertex of the ith triangle ( 0 <= j <= 2 )
				- triang.at(i).at(j).at(k) is the kth coordinate of the jth vertex of the ith triangle ( 0 <= k <= 2 )

So, the parametric modeler class, must have methods `.design_space()`, `.set_design()` and `.triangulate(int N_triangles)` with the behaviour described above. A bare-bones template for this class can be found in `modellers/ParametricModeller.hpp`.

# Description of GMGSA approach

GMGSA [1] is intended for the study of physics-based problems, where the physical Quantities of Interest (QoI) are dependant on the geometric moments of the design under question. Quoting from paragraph 5.10 of [1]:

>Our use of geometric moments is based on the fact that, like most physical quantities, moments are sensitive to the variation of shape features, and the sensitive parameters are those with a high effect on the shape and thus on the associated physics. However, it is not unlikely that some parameters may have a high impact on the shape but a negligible impact on the physics in a design problem ... 
Thus, a good understanding of the underlying physics is necessary to perform a geometric-moment dependent SA

With this in mind, a brief description of the method follows. We assume that we are given some design $\cal{G}$ in $\mathbb{R}^3$, parametrized via $n$ parameters $t_i,\ i=1,...,n$. This means that there is a parametric modeller $\cal{P}$ which for each ${\bf t}\in \mathbb{R}^n$, produces some new shape $\cal{P}({\bf t})$. We can then construct the Shape Signature Vector (SSV) of order $s$ of $\cal{P}({\bf t})$. SSV comprises of all moments of $\cal{P}({\bf t})$ and their invariants up to order $s$. The goal is to measure the effect each parameter has on the SSV, so that a subset of the most important parameters can be identified.

For univariate outputs, Sobol's global sensitivity indices [2] are based on a decomposition of the output-function into functions of progressively more inputs (see eq.12 of [1]) reffered to as ANOVA (functional ANalysis Of VAriance) decomposition. The terms with only one parameter capture the effect of said parameter on the variance of the output. The other terms with more parameters capture the *leftover* effect of the interaction between these parameters on the variance of the output. However, in general, SSV will be an element of $\mathbb{R}^k$, $k>1$ so univariate SA approaches will not suffice. In [3,4] a generalization of this approach to multivariate outputs is developed, called the Covariance-Decomposition Approach.

Finally, having calculated the sensitivity indices of each parameter $t_i$ with respect to the SSV, there are a number of approaches for selecting a subset of parameters, as discussed in [1]. The output of this tool will be the sensitivity indices thus this choice is left up to the user.


# References

1. Shahroz Khan, Panagiotis Kaklis, Andrea Serani, Matteo Diez, _Geometric Moment-Dependent Global Sensitivity Analysis without Simulation Data: Application to Ship Hull Form Optimisation_, 2022

2. Sobol IM. _Global sensitivity indices for nonlinear mathematical models and their Monte Carlo estimates_. Math Comput Simulation 2001;55(1–3):271–80. http://dx.doi.org/10.1016/S0378-4754(00)00270-6.

3. Alexandre Janon, Thierry Klein, Agnès Lagnoux, Maëlle Nodet and Clémentine Prieur, _ASYMPTOTIC NORMALITY AND EFFICIENCY OF TWO SOBOL INDEX ESTIMATORS_

4. Fabrice Gamboa, Alexandre Janon, Thierry Klein, Agn`es Lagnoux, _Sensitivity indices for multivariate outputs_

5. Shahroz Khan, Panagiotis Kaklis, _From regional sensitivity to intra-sensitivity for parametric analysis of free-form shapes: Application to ship design_

6. A new mathematical hull‑form with 10‑shape parameters for evaluation of ship response in waves, Sadaoki Matsui
