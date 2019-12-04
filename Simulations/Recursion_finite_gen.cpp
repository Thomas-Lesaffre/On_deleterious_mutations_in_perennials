#include "header.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <vector>
#include <cmath>
#include <csignal>
#include <algorithm>
#include "mt.h"
#include <cstring>

//#include <gmp.h>
//#include <mpfr.h>

using namespace std;

extern MTRand rnd;

extern FILE * fichierE;
extern FILE * fichierS;

extern FILE * file_mut;
extern FILE * file_car;

// For stopping the program with Ctrl-C:

bool cntl_c_bool = false;
void cntl_c_handler (int bidon)
{
	cntl_c_bool = true;
}

// Parameters :
    // modelv : Specifies which model is to be run, that is, mutations affecting... 1: The maintenance cost; 2: The production cost, 3: Survival.
    // Nv : Population size.
    // Lv : Length of the genetic map (in cM).
    // Sv : Adult survival probability.
    // av : Selfing rate.
	// sv : Phenotypic effect of mutations.
	// maintv : Maintenance cost ($c$ in the main text)
	// prodv : Production cost ($\epsilon$ in the main text)
	// hv : Dominance coefficient of mutations
	// Uv : Mutation rate (per haploid genome).
	// stepv : Frequency at which measurements are written in output files.
	// thrv : Number of generations
	// iteratv : iteration number (in the "parameters.txt" file, the number of iterations to be performed is written)
	// resume_v : Whether or not to look for saved informations. 0: start from zero, 1: look for mutations and characteristics saved in mut_save and car_save.
	
// Recursion :

result recursion(int modelv, int Nv, double Lv, double Sv, double alphav, double sv, double maintv, double prodv, double hv, double Uv, int stepv, double thrv, int iteratv, int resume_v, int save_v, int freq_save_v)
{

    // Naming the output file :

		char nomFichier[256];
		
		stringstream nomF;
		
		if(modelv == 1)
		{
			nomF << "maint" << "_N" << Nv << "_S" << Sv << "_U" << Uv << "_a" << alphav << "_s" << sv << "_maint" << maintv << "_prod" << prodv << "_h" << hv << "_L" << Lv;

			nomF << "_it" << iteratv << ".txt";
			nomF >> nomFichier;
		}
		else if(modelv == 2)
		{
			nomF << "prod" << "_N" << Nv << "_S" << Sv << "_U" << Uv << "_a" << alphav << "_s" << sv << "_maint" << maintv << "_prod" << prodv << "_h" << hv << "_L" << Lv;

			nomF << "_it" << iteratv << ".txt";
			nomF >> nomFichier;		
		}
		else if(modelv == 3)
		{
			nomF << "surv" << "_N" << Nv << "_S" << Sv << "_U" << Uv << "_a" << alphav << "_s" << sv << "_maint" << maintv << "_prod" << prodv << "_h" << hv << "_L" << Lv;

			nomF << "_it" << iteratv << ".txt";
			nomF >> nomFichier;
		}
		else
		{
			nomF << "WRONG_MODEL_SPEC_" << modelv << ".txt";
			nomF >> nomFichier;
		}
		
		int n_save = 0;
		
		// Creating the stream and file names for backup.
		
		char save_mut[256];
		char save_carac[256];
			
		stringstream smut, scarac;
		
		smut << "mut_save";
		scarac << "carac_save";
	
		smut >> save_mut;
		scarac >> save_carac;
		
		ofstream fout(nomFichier); // Stream name for results
		
		fout << "mc"  << " " << "del_rep" << " " << "mw" << endl; // Columns names
		
   	// For stopping prog. with Ctrl-C :
   	
		int stp;
		
	   	cntl_c_bool = false;
		signal(SIGINT, cntl_c_handler);
	
	// Declaring stuff
	
		int twoN = 2*Nv; 
		double twoL = 2*Lv; 
		
		int i, j, k, nb, m, Nmut, step, Ndo, chr1, chr2, ngen;
		double srand, rrand, arand, sum_size, msz, ma, fit, fitj, wHet, wHom, m_count, mean_fit, mc_step, mf_step, ro, rs, dout, dself, delrep_step, wo, ws, Sw, nrand;
		
		vector<double> vec_mc, vec_size, vec_mf, vec_d;
		
		chr * pop = new chr[twoN]; // See what is contained in the "chr" (= "chromosome") structure in the header, "header.h"
		chr * par = new chr[twoN];
		
		chr off1, off2;
		
		result Res;
			
	// Fitness effects 
	
		if(modelv != 3)
		{
			wHet = 1 + hv*sv;
			wHom = 1 + sv;
		}
		else
		{
			wHet = 1 - hv*sv;
			wHom = 1 - sv;		
		}
		
		Ndo = 0;
		step = -1;
		
		vec_mc.clear();
		vec_d.clear();
		
		rs = 0;
		ro = 0;
	
		dself = 0;
		dout = 0;
		
	// Initialization
	
	if(resume_v == 0) // If we start from zero
	{
		for(i = 0; i < Nv; i++)
		{
			nb = 2*i;
		
			pop[nb].age = 1; // Everybody starts at age 1...
			pop[nb+1].age = 1;
		
			pop[nb].sze = growth(pop[nb], pop[nb+1], wHet, wHom, maintv, prodv, modelv); // ... with the same initial size...
			pop[nb+1].sze = growth(pop[nb], pop[nb+1], wHet, wHom, maintv, prodv, modelv);
		
			pop[nb].rep = 0; // ... no offspring produced yet...
			pop[nb+1].rep = 0;
		
			arand = rnd.rand(); // ... depending on the selfing rate, individuals are stochastically set as...
		
			if(arand < alphav)
			{
				pop[nb].S = 1; // ... selfed...
				pop[nb+1].S = 1;
			}
			else
			{
				pop[nb].S = 0; // ... or outcrossed.
				pop[nb+1].S = 0;		
			}
		
			vec_d.push_back(0); // The vector of death indicators is set to zero for everybody
		}
	}
	else // If we start using a backup, we copy informations from backup files to the population vector.
	{
		open_file_mut(); 
		transfer_mut(pop);
	
		open_file_car();
		transfer_car(pop);
		
		for(i = 0; i < Nv; i++)
		{
			vec_d.push_back(0);
		}
	}
	
	// Recursion
		
	for(ngen=0; ngen<thrv; ngen++)
	{	
	
		Ndo++;
		
		// Copying 'pop' to 'par' and building cumulated sizes vector.
		
		vec_size.clear();
		vec_d.clear();
		sum_size = 0;
		
		for(i=0; i<Nv; i++)
		{
			nb = 2*i;
			
			par[nb] = pop[nb];
			par[nb+1] = pop[nb+1];
						
			sum_size += (par[nb].sze + par[nb+1].sze)/2;
			
			vec_size.push_back(sum_size);
			vec_d.push_back(0);
		}

		for(i=0; i<Nv; i++)
		{
			nb = 2*i;
			
			srand = rnd.rand();
			
			Sw = Sv*fitness(pop[nb], pop[nb+1], wHet, wHom, modelv); // Calculate survival probability
			
			if(srand < Sw) // If the individual survives.
			{
				// It grows and ages.
				
				pop[nb].age += 1;
				pop[nb+1].age += 1;
				
				pop[nb].sze = growth(pop[nb], pop[nb+1], wHet, wHom, maintv, prodv, modelv);
				pop[nb+1].sze = growth(pop[nb], pop[nb+1], wHet, wHom, maintv, prodv, modelv);
				
			}
			else // If the individual dies.
			{					
				
				vec_d[i] = 1;
				
				// Build a replacement
				
				do{
				
					// Find a parent
					
						rrand = rnd.rand(sum_size);
						k = -1;
						
						do{
							k++;
						}while(vec_size[k] < rrand);
						
					// Selfing or outcrossing ?
					
						arand = rnd.rand();
						
						if(arand < alphav) // Selfing
						{
							chr1 = 2*k + rnd.randInt(1);
							
							if(chr1 % 2 == 0)
							{
								chr2 = chr1 + rnd.randInt(1);
								par[chr1].rep += 1;
								par[chr1+1].rep += 1;
								
								rec(off1, par[chr1], par[chr1 + 1], Lv);
							
								if(chr2 % 2 == 0)
								{
									rec(off2, par[chr2], par[chr2 + 1], Lv);														
								}
								else
								{
									rec(off2, par[chr2], par[chr2 - 1], Lv);	
								}
							}
							else
							{
								chr2 = chr1 - rnd.randInt(1);
								par[chr1].rep += 1;
								par[chr1-1].rep += 1;
																								
								rec(off1, par[chr1], par[chr1 - 1], Lv);
							
								if(chr2 % 2 == 0)
								{
									rec(off2, par[chr2], par[chr2 + 1], Lv);														
								}
								else
								{
									rec(off2, par[chr2], par[chr2 - 1], Lv);
								}
							}
							
							// Mutation
							
							Nmut = poisdev(Uv);
	
			   				for(m = 0; m < Nmut; m++)
			   				{
			   					off1.sel.push_back(twoL*rnd.rand());
			   					sort(off1.sel.begin(), off1.sel.end());
			   				}

			   				Nmut = poisdev(Uv);
			   				
			   				for(m = 0; m < Nmut; m++)
			   				{
			   					off2.sel.push_back(twoL*rnd.rand());
			   					sort(off2.sel.begin(), off2.sel.end());
			   				}	
			   				
			   				// Setting stuff right (just to be sure !).
							
							off1.age = 1;
							off2.age = 1;
							
							off2.sze = growth(off1, off2, wHet, wHom, maintv, prodv, modelv);
							off1.sze = growth(off1, off2, wHet, wHom, maintv, prodv, modelv);
							
							off1.S = 1;
							off2.S = 1;	
							
							off1.rep = 0;
							off2.rep = 0;
						}
						else // Outcrossing
						{
							chr1 = 2*k + rnd.randInt(1);
							
							if(chr1 % 2 == 0)
							{
								rec(off1, par[chr1], par[chr1 + 1], Lv);
								par[chr1].rep += 1;
								par[chr1+1].rep += 1;							
							}
							else
							{
								rec(off1, par[chr1], par[chr1 - 1], Lv);
								par[chr1].rep += 1;
								par[chr1-1].rep += 1;	
							}
							
							// Choose a father
							
							rrand = rnd.rand(sum_size);
							k = -1;
					
							do{
								k++;
							}while(vec_size[k] < rrand);
							
							chr2 = 2*k + rnd.randInt(1);
								
							if(chr2 % 2 == 0)
							{
								rec(off2, par[chr2], par[chr2 + 1], Lv);
								par[chr2].rep += 1;
								par[chr2+1].rep += 1;
							}
							else
							{
								rec(off2, par[chr2], par[chr2 - 1], Lv);
								par[chr2].rep += 1;
								par[chr2-1].rep += 1;
							}
								
							// Mutation
							
							Nmut = poisdev(Uv);
	
			   				for(m = 0; m < Nmut; m++)
			   				{
			   					off1.sel.push_back(twoL*rnd.rand());
			   					sort(off1.sel.begin(), off1.sel.end());
			   				}

			   				Nmut = poisdev(Uv);
			   				
			   				for(m = 0; m < Nmut; m++)
			   				{
			   					off2.sel.push_back(twoL*rnd.rand());
			   					sort(off2.sel.begin(), off2.sel.end());
			   				}	
								
							// Setting stuff right (just to be sure !)
							
							off1.age = 1;
							off2.age = 1;
							
							off2.sze = growth(off1, off2, wHet, wHom, maintv, prodv, modelv);
							off1.sze = growth(off1, off2, wHet, wHom, maintv, prodv, modelv);
							
							off1.S = 0;
							off2.S = 0;	
							
							off1.rep = 0;
							off2.rep = 0;		
						}
						
						// Incorporate into the population.
						
						nrand = rnd.rand();
						
						if(fitness(off1, off2, wHet, wHom, modelv) < nrand)
						{
							dout += 1 - off1.S;
							dself += off1.S;
						}
						
				}while(fitness(off1, off2, wHet, wHom, modelv) < nrand);
						
				pop[nb] = off1;
				pop[nb+1] = off2;						
			} // End of 'if' 
	
		} // End of 'for' loop over 'pop'
		
		for(i=0; i<Nv; i++)
		{
			nb = 2*i;
			if(vec_d[i] == 0)
			{
				pop[nb].rep = par[nb].rep;
				pop[nb+1].rep = par[nb+1].rep;
			}
			else
			{
				rs += (par[nb].rep*par[nb].S + par[nb+1].rep*par[nb+1].S)/2;
				dself += (par[nb].S + par[nb+1].S)/2;
				
				ro += (par[nb].rep*(1-par[nb].S) + par[nb+1].rep*(1-par[nb+1].S))/2;
				dout += (1 - par[nb].S + 1 - par[nb+1].S)/2;
			}
		}
		
		if (cntl_c_bool) // Control for stopping
		{
			stp = 1;
			break;
		}
		
		m_count = 0;		
		
		mean_fit = 0;
		
		for(i=0; i<Nv; i++)
		{
			nb = 2*i;
			
			m_count += (pop[nb].sel.size() + pop[nb+1].sel.size())/2;
			
			mean_fit += fitness(pop[nb],pop[nb+1],wHet,wHom,3);
		}
		
		m_count /= Nv;
		mean_fit /= Nv;

		vec_mc.push_back(m_count);
		vec_mf.push_back(mean_fit);
		
		// Measures every stepv
		
		if(Ndo % stepv == 0)
		{
			step++;
				
			mc_step = 0;
			delrep_step = 0;
			mf_step = 0;
			
			for(i=0; i<vec_mc.size(); i++)
			{
				mc_step += vec_mc[i];
				mf_step += vec_mf[i];				
			}
			
			mc_step /= vec_mc.size();
			mf_step /= vec_mf.size();
			
			vec_mc.clear();
			vec_mf.clear();
			
			delrep_step = 1 - (rs/dself)/(ro/dout);
			
			rs = 0;
			ro = 0;

			dself = 0;
			dout = 0;
						
			// Appending output file
			
			fout << mc_step << " " << delrep_step << " " << mf_step << endl;
			
			
		}
		
		// Saving population
		
		if((save_v == 1)&&(Ndo % freq_save_v == 0))
		{ 
			if(n_save == 0)
			{
				n_save++;
			}
			else
			{
				remove(save_mut);
				remove(save_carac);
			}
			
			
			ofstream fout_mut(save_mut);
			
			fout_mut.open(save_mut, std::ofstream::out | std::ofstream::trunc);
			fout_mut.close();
			fout_mut.open(save_mut);
			
			ofstream fout_carac(save_carac);
			
			for(i=0; i<twoN; i++)
			{
				
				// fout_mut << "* ";
				// fout_car << "* ";
				
				for(j=0; j<pop[i].sel.size(); j++)
				{
					if(j < pop[i].sel.size() - 1)
					{
						fout_mut << pop[i].sel[j] << " ";
					}
					else
					{
						fout_mut << pop[i].sel[j];
					}
				}
				

				
				if(i < twoN - 1)
				{
					fout_mut << "*";
					fout_carac << pop[i].S << " " << pop[i].age << " " << pop[i].sze << " " << pop[i].rep << "*";
				}
				else
				{
					fout_carac << pop[i].S << " " << pop[i].age << " " << pop[i].sze << " " << pop[i].rep;
				}
			}
			
			cout << "\n Back up completed.\n" << endl;
		}		
		
	}

    delete[] pop;
	delete[] par;
		
    return(Res);
}

