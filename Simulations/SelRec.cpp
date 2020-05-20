// Functions for survival and fitness computation and for recombination:

#include "header.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include "mt.h"
//#include <gmp.h>
//#include <mpfr.h>

using namespace std;

extern MTRand rnd;

 
double fitness(chr &c1, chr &c2, double wHe, double wHo, int model)
{

	double w = 1.0;
	
	if(model == 3)
	{
		int s1 = c1.sel.size() - 1;
		int s2 = c2.sel.size() - 1;

		while ((s1 > -1) || (s2 > -1)) // As long as there is at least one deleterious mutation,
		{
			if (s1 == -1) // If there is no del. mut. on chr1, then all mutations on chr2 are in heterozygous state
			{
				w *= pow(wHe, s2 + 1); //  Thus, fitness is a power function : w = w*(1-hs)^(s2_a+1) where s2_a+1 is nb. of mut. on chr2.

				break;
			}

			if (s2 == -1) // If there is no del. mut. on chr2, same process as before but with the other chr.
			{
				w *= pow(wHe, s1 + 1);
				break;
			}

			// heterozygous mutation on c1:

			if (c1.sel[s1] > c2.sel[s2]) // If the last mutation of chr1 is further along the chr. than that of chr2,
			{
				w *= wHe; // Multiply fitness by (1-hs), that is, add an heterozygous mutation
				s1--; // Diminish s1_a by 1, so that next time, the previous mutation in the vector is analysed.
				continue; // Skip the rest of the loop, go back to top.
			}

			// heterozygous mutation on c2:

			if (c1.sel[s1] < c2.sel[s2]) // Same process but with chr2 this time.
			{
				w *= wHe;
				s2--;
				continue;
			}

			// When we get here, it means we have found c1.B[s1] = c2.B[s2] (none of the 'if' were found true), that is, an homozygous mutation.

			w *= wHo; // Multiply fitness by (1-s)

			s1--; // Diminish mutation count by one on chr1
			s2--; // and on chr2.
		
		} // Back to the top.
	}
	
	return w;
}

double growth(chr &c1, chr &c2, double wHe, double wHo, double c, double e, int model)
{
	
	double g;
	
	if(model != 3)
	{
	
		int s1 = c1.sel.size() - 1;
		int s2 = c2.sel.size() - 1;

		double w = 1.0; // Set fitness to 1.
		
		while ((s1 > -1) || (s2 > -1)) // As long as there is at least one deleterious mutation,
		{
			if (s1 == -1) // If there is no del. mut. on chr1, then all mutations on chr2 are in heterozygous state
			{
				w *= pow(wHe, s2 + 1); //  Thus, fitness is a power function : w = w*(1-hs)^(s2_a+1) where s2_a+1 is nb. of mut. on chr2.

				break;
			}

			if (s2 == -1) // If there is no del. mut. on chr2, same process as before but with the other chr.
			{
				w *= pow(wHe, s1 + 1);
				break;
			}

			// heterozygous mutation on c1:

			if (c1.sel[s1] > c2.sel[s2]) // If the last mutation of chr1 is further along the chr. than that of chr2,
			{
				w *= wHe; // Multiply fitness by (1-hs), that is, add an heterozygous mutation
				s1--; // Diminish s1_a by 1, so that next time, the previous mutation in the vector is analysed.
				continue; // Skip the rest of the loop, go back to top.
			}

			// heterozygous mutation on c2:

			if (c1.sel[s1] < c2.sel[s2]) // Same process but with chr2 this time.
			{
				w *= wHe;
				s2--;
				continue;
			}

			// When we get here, it means we have found c1.B[s1] = c2.B[s2] (none of the 'if' were found true), that is, an homozygous mutation.

			w *= wHo; // Multiply fitness by (1-s)

			s1--; // Diminish mutation count by one on chr1
			s2--; // and on chr2.
		
		} // Back to the top.
		
		if(model == 1)
		{
			long double g1 = 1/(pow(c*w,4.0));
			long double g2 = g1*exp((-(c*w)/e)*((c1.age + c2.age)/2));
			long double g3 = -4*g1*exp((-(3*c*w)/(4*e))*((c1.age + c2.age)/2));
			long double g4 = 6*g1*exp((-(c*w)/(2*e))*((c1.age + c2.age)/2));
			long double g5 = -4*g1*exp((-(c*w)/(4*e))*((c1.age + c2.age)/2));
			g = g1 + g2 + g3 + g4 + g5;
		}
		else if(model == 2)
		{
			long double g1 = 1/(pow(c,4.0));
			long double g2 = g1*exp((-(c)/(e*w))*((c1.age + c2.age)/2));
			long double g3 = -4*g1*exp((-(3*c)/(4*e*w))*((c1.age + c2.age)/2));
			long double g4 = 6*g1*exp((-(c)/(2*e*w))*((c1.age + c2.age)/2));
			long double g5 = -4*g1*exp((-(c)/(4*e*w))*((c1.age + c2.age)/2));
			g = g1 + g2 + g3 + g4 + g5;
		}
		else
		{
			long double g1 = 1/(pow(c,4.0));
			long double g2 = g1*exp((-(c)/e)*((c1.age + c2.age)/2));
			long double g3 = -4*g1*exp((-(3*c)/(4*e))*((c1.age + c2.age)/2));
			long double g4 = 6*g1*exp((-(c)/(2*e))*((c1.age + c2.age)/2));
			long double g5 = -4*g1*exp((-(c)/(4*e))*((c1.age + c2.age)/2));
			g = g1 + g2 + g3 + g4 + g5;
		}	
	}
	else
	{
		long double g1 = 1/(pow(c,4.0));
		long double g2 = g1*exp((-(c)/e)*((c1.age + c2.age)/2));
		long double g3 = -4*g1*exp((-(3*c)/(4*e))*((c1.age + c2.age)/2));
		long double g4 = 6*g1*exp((-(c)/(2*e))*((c1.age + c2.age)/2));
		long double g5 = -4*g1*exp((-(c)/(4*e))*((c1.age + c2.age)/2));
		g = g1 + g2 + g3 + g4 + g5;
	}
	
	return g;
}


// Constructs chromosome "res" by recombining chromosomes c1 and c2
// "Sz" is the map length L (average number of cross-overs per meiosis).
// The S-locus is inherited from chromosome c1.
// The number and position of cross-overs are sampled randomly.

void rec(chr &res, chr &c1, chr &c2, double Sz)
{
	int sz1 = c1.sel.size();
	int sz2 = c2.sel.size();
	
	int s1 = sz1 - 1;
	int s2 = sz2 - 1;
	
	res.sel.clear();
	res.S = 0;
	res.age = 0;
	res.sze = 0;

		int j;
		chr Co;

		double twoS = 2 * Sz;

		// number of cross-overs:

		int nbCo = int(poisdev(Sz));  // Nb. of crossing-overs sampled from a Poisson law with mean "Sz".

		// positions of cross-overs (between 0 and 2L) are put in the vector Co.sel

		for (j = 0; j < nbCo; j++)
        {
            Co.sel.push_back(twoS * rnd.rand()); // Positions of C-O are are sampled at random in twoS.
            sort(Co.sel.begin(), Co.sel.end()); // Positions are sorted.
        }

		// Counting the number of cross-overs on the right of the Modifier locus:

		int cmpt = 0; // set C-O counter to 0.
		j = Co.sel.size() - 1;
		
		while ((j > -1) && (Co.sel[j] > Sz)) // While there is at least 1 C-O and and the considered one is located on right of the M-locus.
		{
			cmpt++; // raise number of C-Os by 1.
			j--; // lower the number of C-O yet to be considered by one.
		}

		// If the number of C-Os on the right is even:

		if (cmpt % 2 == 0)
		{

			for (j = 1; j <= cmpt; j++)
			{
				if (j % 2 == 1) // First C-O results in offspring getting part of chr1, then 2nd C-O (where 'j' is even, see below)
				{               // gives part of chr2, up until the next C-O...
					// All mutations on the right of cross-over on chromosome 1 are incorporated

					while ((s1 > -1) && (c1.sel[s1] > Co.sel[nbCo - j]))
					{
						res.sel.push_back(c1.sel[s1]);
						s1--;
					}
					while ((s2 > -1) && (c2.sel[s2] > Co.sel[nbCo - j]))
						s2--;

				}
				else
				{
					// all mutations on the right of cross-over on chromosome 2 are incorporated

					while ((s2 > -1) && (c2.sel[s2] > Co.sel[nbCo - j]))
					{
						res.sel.push_back(c2.sel[s2]);
						s2--;
					}
					while ((s1 > -1) && (c1.sel[s1] > Co.sel[nbCo - j]))
						s1--;
				}

			}

		}

		// If the number of cross-overs on the right is odd:

		else
		{
			// For each cross-over (starting on extreme right):

			for (j = 1; j <= cmpt; j++)
			{
				if (j % 2 == 1)
				{
					// all mutations on the right of cross-over on chromosome 2 are incorporated

					while ((s2 > -1) && (c2.sel[s2] > Co.sel[nbCo - j]))
					{
						res.sel.push_back(c2.sel[s2]);
						s2--;
					}
					while ((s1 > -1) && (c1.sel[s1] > Co.sel[nbCo - j]))
						s1--;
				}
				else
				{
					// all mutations on the right of cross-over on chromosome 1 are incorporated

					while ((s1 > -1) && (c1.sel[s1] > Co.sel[nbCo - j]))
					{
						res.sel.push_back(c1.sel[s1]);
						s1--;
					}
					while ((s2 > -1) && (c2.sel[s2] > Co.sel[nbCo - j]))
						s2--;
				}

			}

		}


		// mutations between the Modifier locus and the nearest cross-over on the right (on chromosome 1):

		while ((s1 > -1) && (c1.sel[s1] > Sz))
		{
			res.sel.push_back(c1.sel[s1]);
			s1--;
		}
		while ((s2 > -1) && (c2.sel[s2] > Sz))
			s2--;

		// number of cross-overs on the left of the S-locus:

		int frst = nbCo - cmpt;

		// for each cross-over (starting with the nearest to the S-locus on the left):

		for (j = 1; j <= frst; j++)
		{
			if (j % 2 == 1)
			{
				// all mutations on the right of cross-over on chromosome 1 are incorporated

				while ((s1 > -1) && (c1.sel[s1] > Co.sel[frst - j]))
				{
					res.sel.push_back(c1.sel[s1]);
					s1--;
				}
				while ((s2 > -1) && (c2.sel[s2] > Co.sel[frst - j]))
					s2--;
			}
			else
			{
				// all mutations on the right of cross-over on chromosome 2 are incorporated

				while ((s2 > -1) && (c2.sel[s2] > Co.sel[frst - j]))
				{
					res.sel.push_back(c2.sel[s2]);
					s2--;
				}
				while ((s1 > -1) && (c1.sel[s1] > Co.sel[frst - j]))
					s1--;
			}


		}


		// mutations on the left of the left-most cross-over:

		if (frst % 2 == 0)
		{

			while (s1 > -1)
			{
				res.sel.push_back(c1.sel[s1]);
				s1--;
			}

        }
        else
        {

			while (s2 > -1)
			{
				res.sel.push_back(c2.sel[s2]);
				s2--;
			}

        }

	// sorts mutations on offspring chromosome:

	sort(res.sel.begin(), res.sel.end());
}
