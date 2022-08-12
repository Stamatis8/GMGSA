# Geometric Moment-dependent Global Sensitivity Analysis (GMGSA)

## General

- Header-only library

## Description

GMGSA [1] is intended for the study of physics-based problems, where the physical Quantities of Interest (QoI) are dependant on the geometric moments of the design under question. Quoting from paragraph 5.10 of [1]:

>Our use of geometric moments is based on the fact that, like most physical quantities, moments are sensitive to the variation of shape features, and the sensitive parameters are those with a high effect on the shape and thus on the associated physics. However, it is not unlikely that some parameters may have a high impact on the shape but a negligible impact on the physics in a design problem ... 
Thus, a good understanding of the underlying physics is necessary to perform a geometric-moment dependent SA

With this in mind, a brief description of the method follows. We assume that we are given some design $\cal{G}$ in $\mathbb{R}^3$, parametrized via $n$ parameters $t_i,\ i=1,...,n$. This means that there is a parametric modeller $\cal{P}$ which for each ${\bf t}\in \mathbb{R}^n$, produces some new shape $\cal{P}({\bf t})$. We can then construct the Shape Signature Vector (SSV) of order $s$ of $\cal{P}({\bf t})$. SSV comprises of all moments of $\cal{P}({\bf t})$ and their invariants up to order $s$. The goal is to measure the effect each parameter has on the SSV, so that a subset of the most important parameters can be identified.

For univariate outputs, Sobol's global sensitivity indices [2] are based on a decomposition of the output-function into functions of progressively more inputs (see eq.12 of [1]) reffered to as ANOVA (functional ANalysis Of VAriance) decomposition. The terms with only one parameter capture the effect of said parameter on the variance of the output. The other terms with more parameters capture the *leftover* effect of the interaction between these parameters on the variance of the output. However, in general, SSV will be an element of $\mathbb{R}^k$, $k>1$ so univariate SA approaches will not suffice. In [3,4] a generalization of this approach to multivariate outputs is developed, called the Covariance-Decomposition Approach.

Finally, having calculated the sensitivity indices of each parameter $t_i$ with respect to the SSV, there are a number of approaches to selecting a subset of parameters, as discussed in [1]. The output of this tool will be the sensitivity indices thus this choice is left up to the user.

## Implementation



## References

1. Shahroz Khan, Panagiotis Kaklis, Andrea Serani, Matteo Diez, _Geometric Moment-Dependent Global Sensitivity Analysis without Simulation Data: Application to Ship Hull Form Optimisation_, 2022

2. Sobol IM. _Global sensitivity indices for nonlinear mathematical models and their Monte Carlo estimates_. Math Comput Simulation 2001;55(1–3):271–80. http://dx.doi.org/10.1016/S0378-4754(00)00270-6.

3. Alexandre Janon, Thierry Klein, Agnès Lagnoux, Maëlle Nodet and Clémentine Prieur, _ASYMPTOTIC NORMALITY AND EFFICIENCY OF TWO SOBOL INDEX ESTIMATORS_

4. Fabrice Gamboa, Alexandre Janon, Thierry Klein, Agn`es Lagnoux, _Sensitivity indices for multivariate outputs_
