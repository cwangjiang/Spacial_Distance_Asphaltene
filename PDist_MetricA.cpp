#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <fstream>
#include <cmath>

using namespace std;

const int N = 8;//total # of molecules
const int allB = 21;//total # of bead for each molecule

int main()
{
int T = 11;//total snapshot want, /ns
int NB = 21;//number of bead want to include
int exNB = 0;//number of bead not included 
//int N = 200;//molecule number
double A[N][allB][3] = {0.0};//coordinate matrix updated at each time
double r[3]={0.0};//temperal bead to bead displancement
double r2 = 0.0;//temperal bead to bead distance
double PD = 0.0;//molecule to molecule distance is the minima of al bead to bead distance
double L = 7.2;//box size, we need to change it for difference simulation
double min = 7.2;//initial minima bead to bead distance
string line; //temperal string for remove extra line in .gro file


FILE* Distance = fopen("Distance_MetricA.dat","w");//open pairwise distance matrix file

ifstream in("mappedALL.gro");//open original .gro coordinate data file


if(!in)
{
cout<<"Fail of open file"<<endl;
return 0;
}

for(int t = 0 ; t <T ; t++){
	cout<<t<<endl;//output frame time to see processes
	getline(in,line);//remove first two lines (frame time and # of bead) for each frame
	getline(in,line);
	for(int m = 0 ; m < N ; m++){
		for(int b = 0 ;b <allB ; b++){
			in.ignore(20);//skip first colums for each line
			for(int c = 0 ; c<3 ; c++){
				in>>A[m][b][c];//fill coordinate matrix for each frame
			}
			getline(in,line);//remove final part for each line, so that we can retrun to next line
		}
	}
	getline(in,line);//remove the last line (box size) for each frame

//conpute pariwise distance
for(int m1 = 0; m1 < N ; m1++){//loop over m1 m2 all combination of molecules
	for(int m2 = m1+1; m2 < N ; m2++){//we just need to compute halp number of melecule molecule pair, since melecule molecule distance are symmetrical
		for(int i = 0; i < NB ; i++){//loop over all possible combination of bead to bead pair in such molecule pair
			for(int j = 0; j < NB ; j++){
				for(int k = 0; k < 3 ; k++){//standard way of computeing distance in periodic box
					r[k] = A[m1][i][k] - A[m2][j][k];
					r[k] = r[k] - L*round(r[k]/L);
					r2 += r[k]*r[k];
				}
				r2 = sqrt(r2);//distace between beads
				if (r2 < min){
				min = r2;//find the minima bead to bead in a molecular pair
				}
			r2 = 0.0;//reset bead to bead distance and move to next bead to bead pair
			}
		}
	PD = min;//after loop overn all bead pair in a molecular pair, we got the molecule to molecule distance.
	fprintf(Distance,"%f\n",PD);//output molecular molecular pairwise distance to file
	min = L;//reset min to a large value
	}
}
}
}

