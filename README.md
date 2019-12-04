# On Deleterious mutations in perennials


This repository relates to a manuscript that will soon be posted on _BioRXiv_. 

In this paper, we ask whether mutations affecting fitness differently with respect to life-history may explain the large differences in inbreeding depression observed between perennials and annuals. To do so, we combine a physiological growth model and multilocus population genetics approaches in order to describe a full genotype-to-phenotype-to-fitness map, where the phenotype relates to fitness through biological assumptions, so that the fitness landscape emerges from biological assumptions instead of being assumed a priori. We study the behaviour of different types of mutations affecting growth or survival, and explore their consequences in terms of inbreeding depression and mutation load. Then, we discuss the role deleterious mutations maintained at mutation-selection balance may play in the coevolution between growth and survival strategies.

The _Simulations_ folder contains our individual-centered simulations code. The _parameters.txt_ file is a template for the parameters input file. The files necessary for program compiling are listed below.
* _fichiers.cpp_ contains the functions relating to files management in general (reading parameters, managing backups...).
* _main.cpp_ is the main file of the program. It calls all the other functions used.
* _ranbin.cpp_ contains the functions we use to simulate various distributions (mostly Poisson) using a random numbers generator.
* _SelRec.cpp_ contains the recombination, growth and fitness calculation functions.
* _Recursion_finite_gen.cpp_ is the most important file of this program. It contains the "recursion()" function, which performs the simulations. 
* _mt.h_ is the MersenneTwister random numbers generator.
* _header.h_ contains all the prototypes of functions and structures used in the program.

The _Numerical calculations_ folder contains the _Mathematica_ notebooks we used for numerical calculations. Because we use two different methods to obtain numerical results, there are two notebooks, _Approximations_lf.nb_ for the Lifetime Fitness approach, and _Approximations_lc.nb_ for the Life Cycle approach.

If readers need any explanation regarding the programs and notebooks presented in this repository, they can e-mail me at: *thomas.lesaffre@univ-lille.fr*
