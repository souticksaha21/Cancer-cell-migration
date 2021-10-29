//Code for multiple reactions on a circle A\toA*,A*\to A, using gillespie Algo
#include "stdio.h"
#include "math.h"

double dt=5;//0.0005;
double T=36; //total time

int NT=(int)(T/dt);
double k1=0.1;
int MAXTRIAL=400; //number of trials or trajectories
int itbunch=100; //how many dt's you want to bunch in the average
int itbmax=(int)(NT/itbunch)+1;
double pi=3.14159;
int dNtheta=50; //number of divisions of the circle
double dtheta=2*pi/(float)(dNtheta);
double a=20;//size of the cell in um
double R=0.5*a;//radius of the cell
int Amax=100000; //total number of binding molecules
double dr=1;//determine the dr from diffusivity later
double TGFE=25*0.6;//per um^3 conc of TGF
double EGFE=0*0.6;//per um^3 conc of EGF
double gT=50*0.6/1000;//per um^4 gradient of TGF
double gE=0*0.6/1000;//per um^4 gradient of EGF
double eta=0.5;
double beta=0.1;
double V=R*dtheta*pow(dr,2);//relevant volume
double EGF=400*0.6;// for definition of kE
double TGF=25*0.6;// for definition of kT
double kE=1/(beta*V*EGF);
double kT=eta*kE*EGF/TGF;
double velocity=20;//speed in um/s
double persistence=0.0; //persistence in movement 0 no persistence 1 full persistence
int printtraj=10; //trajectory printing per printtraj trials
double eq_time=0.0*T;//equilibriation time
double backgound=1; //if background 1 then the background for each cell is variable, if it is 0 then all cells have the same background EGFE/TGFE.
double motion_effect=1; //if motion_effect 1 then the background for each cell is changes as they move, if it is 0 then background for each cell is does not change as they move
double dist_max=450;//maximum length of length of the channel where cells can lie. cells lie between -dist_max to dist_max
