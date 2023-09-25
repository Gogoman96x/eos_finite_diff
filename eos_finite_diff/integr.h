//Header um Integrale numerisch zu berechnen

#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>

using namespace std;

/*der funktion werde die Parameter für das Integral, die Funktionsparameter und die 
Funktion selbst übergeben.*/

/*Mittels der Stützstellen wird das Integral berechnet mit einem simplen Algo.. Dieser
  Algo. wird noch durch ein besseren ersetzt.*/

double integr(double x1,double x2, double x3, double x, double y,
              double z, double (*func)(double,double,double,double)){//x1=xmin,x2=xmax,x3=stepsx

    double I=0;//Wert des Integrals
    double xn=x1;//Laufvariable für den Integrationsparameter
    double delx=(x2-x1)/x3;                                                                           
    double yn=func(xn,x,y,z); 


    I+=1/2.*delx*yn;
    

    for(int i=0;i<x3;i++){
        xn+=delx;
        yn=func(xn,x,y,z);

        I+=delx*yn;
        
        if(i==x3-1){
            I+=1/2.*delx*yn;
            
            break;
        }
    }
    return I;
}

