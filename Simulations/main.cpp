#include "header.h"
#include "mt.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
//#include <gmp.h>
//#include <mpfr.h>


using namespace std;

// Random number generator:

MTRand rnd;

// Pointers on input and output files:

FILE * fichierE;
FILE * fichierS;

// Pointers on back ups

FILE * file_mut;
FILE * file_car;

int main()
{

    cout << "Program initialization" << "\n";

	// Parameters:

	int N, it, iteration, step, i, resume, save, freq_save, model;
	double L, S, alpha, s, maint, prod, h, U, thr;

	result Rslt;

	// Opens input and output files:

	bool end;
	ouvrirFichierE();
	ouvrirFichierS();
	end = false;

	int no = 1;

	do
	{
        //reads parameter values from input file:

		end = lireFichier(model, N, L, S, alpha, s, maint, prod, h, U, step, thr, it, resume, save, freq_save);
					   	 // end = true if end of input file
        if(!end) cout << "\n___________________________________\n" << "\n";

        if(!end)
            for (i = 0; i < it; i++)
            {

            iteration = i;

            cout << "\n" << "Iteration Number : "<< i << "\n";

			// Simulation:

 			Rslt = recursion(model, N, L, S, alpha, s, maint, prod, h, U, step, thr, iteration, resume, save, freq_save);
 			
 			stringstream command;
			command << "cp mut_save mut_model" << model << "_N" << N << "_S" << S << "_U" << U << "_a" << alpha << "_s" << s << "_maint" << maint << "_prod" << prod << "_h" << h << "_L" << L << "_it" << iteration << ";" << "cp carac_save car_model" << model << "_N" << N << "_S" << S << "_U" << U << "_a" << alpha << "_s" << s << "_maint" << maint << "_prod" << prod << "_h" << h << "_L" << L << "_it" << iteration;
			
			system(command.str().c_str()); 
			
			command.str("");
			
            }

            no++;

	}while(!end);

	// Closes files:
	
    fclose(fichierE);
	fclose(fichierS);
	
	// Tidy up repository
	 
	system("mkdir back_up; mv mut_model* ./back_up/; mv car_model* ./back_up/; rm results.txt");
	 
	if(resume == 1)
	{
		fclose(file_mut);
		fclose(file_car);
	}
	
	return 0;
}


