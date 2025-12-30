#include <math.h>
#include <iostream>
#include "problemi.h"
using namespace std;
typedef double Real;
const Real pi=4.0*atan(1.0);
// ********************************************************
void f_progetto3(Real F[],Real t,Real U[])
{
    F[0]=U[1]*U[2]*sin(t)-U[0]*U[1]*U[2];
    F[1]=-U[0]*U[2]*sin(t)+1.0/20.0*U[0]*U[2];
    F[2]=U[0]*U[0]*U[1]-1.0/20.0*U[0]*U[1];

}

void dati_progetto3(Real *t0,Real *T,Real v[])
{
    *t0=0.0;
    *T=1.0;  //Per disegnare con Gauss3 la traiettoria passante per il punto v sostituisco T=1000
             //(così visualizzo sulla sfera un intervallo di tempo maggiore)

    Real a=1.0/sqrt(3);
    v[0]=a;
    v[1]=a;
    v[2]=a;
}

void dati_progetto3s(Real *t0,Real *T)
{
    *t0=0.0;
    *T=100.0;
}

void Jf_progetto3(Real J[], Real t,Real *U)
{
    J[0]=-U[1]*U[2];
    J[1]=U[2]*sin(t)-U[0]*U[2];
    J[2]=U[1]*sin(t)-U[0]*U[1];
    J[3]=-U[2]*sin(t)+1.0/20.0*U[2];
    J[4]=0;
    J[5]=-U[0]*sin(t)+1.0/20.0*U[0];
    J[6]=2*U[0]*U[1]-1.0/20.0*U[1];
    J[7]=U[0]*U[0]-1.0/20.0*U[0];
    J[8]=0;
}
