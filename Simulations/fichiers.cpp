// Functions to open input and output files,
// read parameter values and write them in output file

#include "header.h"
#include <iostream>
#include <fstream>
#include <sstream>

//#include <gmp.h>
//#include <mpfr.h>

using namespace std;

extern FILE * fichierE;
extern FILE * fichierS;

extern FILE * file_mut;
extern FILE * file_car;

//Opens input file:

void ouvrirFichierE()
{
	fichierE = fopen(fichierLecture,"r");
	if (!fichierE)
		cout << "The file " << fichierLecture
			 << " doesn't exist!\n\n";
}

//Opens output file:

void ouvrirFichierS()
{
	fichierS = fopen(fichierEcriture,"a");
	if (!fichierS)
		cout << "Impossible to open " << fichierEcriture << " !\n\n";
}



// Reads parameter values from input file,
// returns 1 if it reaches the end of input file, else returns 0:

bool lireFichier(int &modr, int &Nr, double &Lr, double &Sr, double &alphar, double &sr, double &maintr, double &prodr, double &hr, double &Ur, int &stepr, double &thrr, int &iterr, int &resume_r, int &save_r, int &freq_save_r)
{
	int z;
	bool term;
	do {z = fgetc(fichierE);} while (!((z == '*') || (z == EOF)));
		// Lines with parameter sets must begin with *
	if (z == EOF)
	{
		cout << "\nEnd of input file\n";
		term = true;
	}
	else
	{
		fscanf(fichierE," %d",&modr);		
		fscanf(fichierE," %d",&Nr);
		fscanf(fichierE," %lf",&Lr);
		fscanf(fichierE," %lf",&Sr);
		fscanf(fichierE," %lf",&alphar);
		fscanf(fichierE," %lf",&sr);
		fscanf(fichierE," %lf",&maintr);
		fscanf(fichierE," %lf",&prodr);
		fscanf(fichierE," %lf",&hr);
		fscanf(fichierE," %lf",&Ur);
		fscanf(fichierE," %d",&stepr);				
		fscanf(fichierE," %lf",&thrr);						
		fscanf(fichierE," %d",&iterr);
		fscanf(fichierE," %d",&resume_r);								
		fscanf(fichierE," %d",&save_r);	
		fscanf(fichierE," %d",&freq_save_r);																					
        term = false;
	}
	
	return term;
}

// Open back up files

void open_file_mut()
{
	file_mut = fopen(fileMut,"r");
	
	if(!file_mut)
	{
		cout << "The file " << file_mut << " doesn't exist !" << endl;
	}
}

void open_file_car()
{
	file_car = fopen(fileCar,"r");
	
	if(!file_car)
	{
		cout << "The file " << file_car << " doesn't exist !" << endl;
	}
}

// Read back up files for mutations

void transfer_mut(chr * v)
{	
	int z, c, k;
	double m;
	char ch;
	
	stringstream tmp;
	
	c = 0;
	
	do{
		tmp.clear();	
		m=0;

		z = fgetc(file_mut);
		
		if(z != EOF)
		{
				ch = z;
		
			if((ch != ' ')&&(ch != '*')&&(z != EOF))
			{
				tmp << ch;
			}
			else if((ch == ' ')&&(z != EOF))
			{
				tmp >> m;
				tmp.str("");
			
				v[c].sel.push_back(m);
			}
			else if((ch == '*')&&(z != EOF))
			{
				tmp >> m;
				tmp.str("");
			
				v[c].sel.push_back(m);
				
//				for(k=0; k<v[c].sel.size(); k++)
//				{
//					cout << v[c].sel[k] << " ";
//				}
//				
//				cout << endl;
				
				c++;
			}
		}
		else
		{
			tmp >> m;
			tmp.str("");
			
			v[c].sel.push_back(m);
		}	
			
	}while(z != EOF);
	
	cout << "Mutations transfered." << endl;
}

// Read back up files for characteristics

void transfer_car(chr * v)
{	
	int z, c, k, car_c;
	double value;
	char ch;
	
	stringstream tmp;
	
	c = 0;
	car_c = 0;
	
	do{
		tmp.clear();	
		value=0;

		z = fgetc(file_car);
		
		if(z != EOF)
		{
			ch = z;
		
			if((ch != ' ')&&(ch != '*')&&(z != EOF))
			{
				tmp << ch;
			}
			else if((ch == ' ')&&(z != EOF))
			{
				tmp >> value;
				tmp.str("");
				
				if(car_c == 0)
				{
					v[c].S = value;
				}
				else if(car_c == 1)
				{
					v[c].age = value;
				}
				else if(car_c == 2)
				{
					v[c].sze = value;
				}
				else if(car_c == 3)
				{
					v[c].rep = value;
				}
				
				car_c++;
			}
			else if((ch == '*')&&(z != EOF))
			{
				tmp >> value;
				tmp.str("");
			
				v[c].S = value;
				
//				for(k=0; k<v[c].sel.size(); k++)
//				{
//					cout << v[c].sel[k] << " ";
//				}
//				
//				cout << endl;
				
				c++;
				car_c = 0;
			}
		}
		else
		{
			tmp >> value;
			tmp.str("");
			
			v[c].sel.push_back(value);
		}	
			
	}while(z != EOF);
	
	cout << "Characteristics transfered." << endl;
}
