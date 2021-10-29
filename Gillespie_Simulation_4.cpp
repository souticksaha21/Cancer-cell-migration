//Code for multiple reactions on a circle A\toA*,A*\to A, using gillespie Algo
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "mersenne_twister.hpp"
#include <time.h>
#include "Gillespie_Simulation_4_param.cpp" //parameter file

int main(void)
{
    time_t timer; //calculating time from local clock for random number generator
    struct tm y2k = {0};
    double seconds; //seed for random number generator from local clock
    
    y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
    y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;
    
    time(&timer);  /* get current time; same as: timer = time(NULL)  */
    
    seconds = difftime(timer,mktime(&y2k));
    MTRand rg(seconds); //initiate random number generator from clock
    FILE *fpA,*fpAvg,*fpcos;
    char filename[FILENAME_MAX];
    double r1,r2;
    int it; //time index
    double A;
    int trials;
    int nit[itbmax];
    int Aavg[itbmax];
    int itb,j;
    double time;
    double tau,p2;
    double alpha0; //total propensity
    int ntheta[Amax];
    double Amol[dNtheta];//number of A molecules in a given dtheta interval
    double prop[dNtheta];//propensity index at dtheta interval
    double prob[2*dNtheta];//probability of reactions
    double probtot[2*dNtheta];//total probability of reactions
    double probst[dNtheta];//probability of Amolst
    double probtotst[dNtheta];//total probability of Amolst
    double Amolst[dNtheta];//number of A* molecules in a given dtheta interval
    double xst,yst;
    double thetast[Amax];
    double Amolsttot;
    double meanCI;
    double thetaold;
    double T0,E0;
    double thetamov; //movement direction of cell via maximum A*
    double dx;
    sprintf(filename,"Gillespie_Simulation_4_T_A0_100000.dat");
    fpA = fopen(filename,"w");
    sprintf(filename,"Gillespie_Simulation_4_CI_T_A0_100000.dat");
    fpAvg = fopen(filename,"w");
    sprintf(filename,"Gillespie_Simulation_4_cos_T_A0_100000.dat");
    fpcos = fopen(filename,"w");


    meanCI=0.0;
    for (trials=1; trials<=MAXTRIAL; trials++) {
            thetaold=2*pi*rg.randDblExc();
        for (int k1=1; k1<=dNtheta; k1++) {
            Amol[k1]=0;
            Amolst[k1]=0;
        }
        for (int j=1; j<=Amax; j++) {//initialising the position of Amax number of A molecules
            double theta=2*pi*rg.randDblExc();//randomly choosing a theta for A molecule
            ntheta[j]=(int)(theta/dtheta)+1;
            int k=(int)(theta/dtheta)+1;
            Amol[k]=Amol[k]+1; //number of A molecules in a given dtheta element
            thetast[j]=0;
        }
        
        if (backgound==1) {
            double t11=rg.randDblExc();
            T0=TGFE+gT*dist_max*(2*t11-1); //starting TGF conc
            E0=EGFE+gE*dist_max*(2*t11-1); //starting EGF conc
        } else {
            T0=TGFE; //starting TGF conc
            E0=EGFE; //starting EGF conc
        }

        
        time=0;
        xst=0;
        yst=0;
        Amolsttot=0.0;
        while (time<T) {
            for (int k1=1; k1<=dNtheta; k1++) {
                prop[k1]=(kT*(T0+gT*R*cos(k1*dtheta))+kE*(E0+gE*R*cos(k1*dtheta)))*V*Amol[k1];
                //printf("%f %f \n",time, Amolst[k1]/prop[k1]);
            }
            alpha0=0;
            for (int k2=1; k2<=dNtheta; k2++) {
                alpha0=alpha0+prop[k2]+Amolst[k2];
                prob[2*k2-1]=prop[k2];
                prob[2*k2]=Amolst[k2];
            }
            
            for (int k31=1; k31<=2*dNtheta; k31++) {
                probtot[k31]=0;
            }
            
            for (int k3=1; k3<=2*dNtheta; k3++) {
                prob[k3]=prob[k3]/alpha0;
                if (k3==1) {
                    probtot[k3]=prob[k3];
                } else {
                    probtot[k3]=probtot[k3-1]+prob[k3];
                }
    
            }
            
            r1=rg.randDblExc();
            r2=rg.randDblExc();
            tau=(log(1/r1))/(alpha0);
            if (r2>0 && r2<=probtot[1]) {
                if (Amol[1]>0) {
                    Amol[1]=Amol[1]-1;
                    Amolst[1]=Amolst[1]+1;
                    Amolsttot=Amolsttot+1;
                    //double r20=rg.randDblExc();
                    //if (r20<persistence && time>0) {
                      //  xst=xst+velocity*tau*cos(thetaold);
                       // yst=yst+velocity*tau*sin(thetaold);
                    //} else {
                      //  xst=xst+velocity*tau*cos(dtheta);
                       // yst=yst+velocity*tau*sin(dtheta);
                        //thetaold=dtheta;
                    //}

                }
                } else {
                for (int k4=1; k4<=2*dNtheta-1; k4++) {
                    if (r2>probtot[k4] && r2<=probtot[k4+1]){
                        if (((k4+1) % 2) == 0) {
                            int k41=(int)(0.5*(k4+1));
                            if (Amolst[k41]>0) {
                                Amolst[k41]=Amolst[k41]-1;
                                Amolsttot=Amolsttot-1;
                                Amol[k41]=Amol[k41]+1;
                                double dth=2*pi*rg.randDblExc();
                                //double r21=rg.randDblExc();
                                //if (r21<persistence && time>0) {
                                  //  xst=xst+velocity*tau*cos(thetaold);
                                    //yst=yst+velocity*tau*sin(thetaold);
                                //} else {
                                  //  xst=xst+velocity*tau*cos(dth);
                                    //yst=yst+velocity*tau*sin(dth);
                                    //thetaold=dth;
                                //}

                            }
                        } else {
                            int k42=(int)(0.5*k4+1);
                            if (Amol[k42]>0) {
                                Amolst[k42]=Amolst[k42]+1;
                                Amolsttot=Amolsttot+1;
                                Amol[k42]=Amol[k42]-1;
                                //double r22=rg.randDblExc();
                                //if (r22<persistence && time>0) {
                                  //  xst=xst+velocity*tau*cos(thetaold);
                                  //  yst=yst+velocity*tau*sin(thetaold);
                                //} else {
                                  //  xst=xst+velocity*tau*cos(k42*dtheta);
                                  //  yst=yst+velocity*tau*sin(k42*dtheta);
                                  //  thetaold=k42*dtheta;
                                //}


                            }
                        }
                    }
                    }
                }
            
            
            for (int k25=1; k25<=dNtheta; k25++) {
                probtotst[k25]=0;
                if (Amolsttot>0) {
                    probst[k25]=Amolst[k25]/Amolsttot;
                } else {
                    probst[k25]=0;
                }
                
                if (k25==1) {
                    probtotst[1]=probst[k25];
                } else {
                    probtotst[k25]=probtotst[k25-1]+probst[k25];
                }
                //printf("k is %d %f \n",k25,probtotst[k25]);
                if (trials==20  && ((int)(time/dt) % 1000) == 0) {
                    fprintf(fpcos,"%f %f \n",cos(k25*dtheta),probst[k25]);
                }
            }
            
            double probmax=0;
            for (int k22=1; k22<=dNtheta; k22++) {
                if (probst[k22]>probmax) {
                    thetamov=k22*dtheta;
                    probmax=probst[k22];
                }
            }
            
            //if (time>eq_time) {
              //  xst=xst+velocity*tau*cos(thetamov);
              //  yst=yst+velocity*tau*sin(thetamov);
            //}
            
            
            double r3=rg.randDblExc();
            if (r3>0 && r3<=probtotst[1]) {
                double r22=rg.randDblExc();
                if (r22<persistence && time>0) {
                    xst=xst+velocity*tau*cos(thetaold);
                    yst=yst+velocity*tau*sin(thetaold);
                    dx=velocity*tau*cos(thetaold);
                } else {
                    xst=xst+velocity*tau*cos(dtheta);
                    yst=yst+velocity*tau*sin(dtheta);
                    //printf("sin1 of angle is %f \n",sin(dtheta));
                    thetaold=dtheta;
                    dx=velocity*tau*cos(dtheta);
                }
            } else {
                for (int k4=1; k4<=dNtheta-1; k4++) {
                    if (r3>probtotst[k4] && r3<=probtotst[k4+1]){
                        double r22=rg.randDblExc();
                        if (r22<persistence && time>0) {
                            xst=xst+velocity*tau*cos(thetaold);
                            yst=yst+velocity*tau*sin(thetaold);
                            dx=velocity*tau*cos(thetaold);
                        } else {
                            xst=xst+velocity*tau*cos((k4+1)*dtheta);
                            yst=yst+velocity*tau*sin((k4+1)*dtheta);
                            thetaold=(k4+1)*dtheta;
                            dx=velocity*tau*cos((k4+1)*dtheta);
                            //printf("sin2 of angle is %f \n",sin((k4+1)*dtheta));
                        }
                    }
                }
            }
            
            if (motion_effect==1) {
                T0=T0+gT*dx;
                E0=E0+gE*dx;
            }

            


            time=time+tau;
            //printf("%f %f \n",time,alpha0);
            if ((trials % printtraj)==0 && ((int)(time/dt) % 10) == 0) {
               fprintf(fpA,"%d %f %f \n",trials,xst,yst);
            }


        }
        if ((trials % printtraj)==0) {
        fprintf(fpA,"\n");
        }
    double dist=pow(xst,2)+pow(yst,2);
    double CI=xst/(pow(dist,0.5));
    printf("%d %f %f %f \n",trials,CI,T0,E0);
        meanCI=meanCI+CI/((float)MAXTRIAL);


    fprintf(fpAvg,"%d %f %f %f %f \n",trials,rg.randDblExc(),CI,T0,E0);
    }
        printf("Mean of CI is %f \n",meanCI);
    fclose(fpA);
    fclose(fpAvg);
    fclose(fpcos);
    return 0;
}
