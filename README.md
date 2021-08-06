# On Deleterious mutations in perennials

This repository relates to the article "On deleterious mutations in perennials: inbreeding depression, mutation load and life-history evolution" published in _The American Naturalist_ (https://doi.org/10.1086/713499).

Zenodo DOI: [![DOI](https://zenodo.org/badge/225906219.svg)](https://zenodo.org/badge/latestdoi/225906219)


In this paper, we ask whether mutations affecting fitness differently with respect to life-history may explain the large differences in inbreeding depression observed between perennials and annuals. To do so, we combine a physiological growth model and multilocus population genetics approaches in order to describe a full genotype-to-phenotype-to-fitness map, where the phenotype relates to fitness through biological assumptions, so that the fitness landscape emerges from biological assumptions instead of being assumed a priori. We study the behaviour of different types of mutations affecting growth or survival, and explore their consequences in terms of inbreeding depression and mutation load. Then, we discuss the role deleterious mutations maintained at mutation-selection balance may play in the coevolution between growth and survival strategies.

The _Simulations_ folder contains our individual-centered simulations code. The _parameters.txt_ file is a template for the parameters input file. The files necessary for program compiling are listed below.
* _fichiers.cpp_ contains the functions relating to files management in general (reading parameters, managing backups...).
* _main.cpp_ is the main file of the program. It calls all the other functions used.
* _ranbin.cpp_ contains the functions we use to simulate various distributions (mostly Poisson) using a random numbers generator.
* _SelRec.cpp_ contains the recombination, growth and fitness calculation functions.
* _Recursion_finite_gen.cpp_ is the most important file of this program. It contains the "recursion()" function, which performs the simulations. 
* _mt.h_ is the MersenneTwister random numbers generator.
* _header.h_ contains all the prototypes of functions and structures used in the program.

The _Numerical_calculations_ folder contains the _Wolfram Language_ scripts we used for numerical calculations. Because we use two different methods to obtain numerical results, there are two scripts, _Approximations_lf.wls_ for the Lifetime Fitness approach, and _Approximations_lc.wls_ for the Life Cycle approach. The files _math_par_lf.txt_ and _math_par_lc.txt_ are the corresponding parameter files, and the _Approximations_lf.nb_ and _Approximations_lc.nb_ are the _Mathematica_ notebook versions of the scripts, which are easier to read because of the formatting. The file _Approximations_lf.nb_ is annotated.

If readers need any explanation regarding the programs and notebooks presented in this repository, they can e-mail me at: *thomas.lesaffre@univ-lille.fr*
