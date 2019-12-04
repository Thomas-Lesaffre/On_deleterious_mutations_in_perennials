#ifndef DEPRESSION_H
#define DEPRESSION_H

#include <vector>
#include <iostream>
// #include "MersenneTwister.h"
#include "mt.h"
//#include <gmp.h>
//#include <mpfr.h>

using namespace std;

// global variables

#define fichierLecture "parameters.txt"     // names of input
#define fichierEcriture "results.txt"		// and output files

#define fileMut "mut_save" // back up files.
#define fileCar "carac_save" 

// definition of structure "chr" representing a chromosome:
// "M" is the Modifier locus
// "sel" is a vector containing the positions of deleterious alleles along the chromosome

struct chr
{
      vector<double> sel; // Selected loci
      
      double S; // 1 = Selfed, 0 = Outcrossed
      double age; // Age of the individual
      double sze; // Size of the individual
      double rep; // Nb. of offspring produced by the individual.

      chr(){}; // Default constructor
      ~chr(){}; // Default destructor
};

struct result
{
       double pouet;
};


// Prototypes of functions

void ouvrirFichierE();

void open_file_mut();
void open_file_car();

void transfer_mut(chr * v);
void transfer_car(chr * v);

void ouvrirFichierS();

bool lireFichier(int &modelr, int &Nr, double &Lr, double &Sr, double &alphar, double &sr, double &maintr, double &prodr, double &hr, double &Ur, int &stepr, double &thrr, int &iterr, int &resume_r, int &save_r, int &freq_save_r);

result recursion(int modelv, int Nv, double Lv, double Sv, double alphav, double sv, double maintv, double prodv, double hv, double Uv, int stepv, double thrv, int iteratv, int resume_v, int save_v, int freq_save_v);

double gammln(const double xx);
double poisdev(const double xm);
double binldev(const double pp, const int n);
double fitness(chr &c1, chr &c2, double wHe, double wHo,  int model);

double growth(chr &c1, chr &c2, double wHe, double wHo, double maint, double prod, int model);

void rec(chr &res, chr &c1, chr &c2, double sz);
void cntl_c_handler(int bidon);


#endif
