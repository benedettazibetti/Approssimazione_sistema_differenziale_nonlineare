// PROGETTO ESAME N3
// METODO: MULTISTEP ESPLICITO
// Utilizzo ADAMS BASHFORTH 4


#include <iostream>
#include <math.h>
#include <fstream>
#include "..\funzioni\problemi.h"
#include "..\funzioni\mat_vet.h"
#include "..\funzioni\Runge_Kutta.h"
#include "..\funzioni\stampe.h"
#include "..\funzioni\multi_step_metodi.h"

using namespace std;
typedef double Real;
const Real pi=4.0*atan(1.0);

// definisco puntatori al tipo di funzioni che enunciano il problema
void(*effe)(Real*, Real,Real*);
void(*dati)(Real *,Real *,Real *);
void(*butcher)(int, Real *,Real *,Real *,Real *, int *); //Carica i coefficienti della matrice di Butcher
void(*Jf)(Real*, Real,Real*);
void (*set_multi_step)(int ,Real*,Real *);

// variabili comuni
Real t,T,h;

// dati del problema, del metodo e dimensioni
const int d=3;
Real u[d];

// dimensiono matrice e vettori Butcher
    const int ns=6;  //Numero degli stadi (ns=6 perchè considero FEHL5)
    const int passi=5;
    Real sigma[passi+1];
    Real ro[passi+1];
    Real Old_U[passi][d];
    Real Old_f[passi][d];
    Real b[ns];
    Real c[ns];
    Real A[ns][ns];
    Real K[ns][d];
    Real Z[d];
    int valf=0;

int main ()
{
    effe=f_progetto3;
    dati=dati_progetto3;
    set_multi_step = AB;
    butcher=FEHL5;
    butcher(ns,A[0],b,c,0,0);

    int ordine=passi;

    // predispongo file stampa
    char n_file[21]= {0};
    cout << "dammi nome file di stampa(max 20 caratteri)";
    cin >> n_file ;
    ofstream prt(n_file);
    prt.precision(14);



    // input
    dati(&t,&T,u);
    set_multi_step(ordine,ro,sigma);
    // Jf=0 se metodo è esplicito
    unsigned long N; //intero privo di segno
    cout << "inserisci numero di passi N  = ";
    cin >> N;
    h=(T-t)/Real(N);

    stampa(d,t,u,&prt);
    stampa(d,t,u);


    // Calcolo dei valori iniziali aggiuntivi con un metodo RK esplicito (ciclo ridotto sul tempo)
    for(unsigned long n=0; n<(passi-1); n++)
    {
        copia(Old_U[n],u,d);
        effe(Old_f[n],t,u);
        copia(K[0],Old_f[n],d);

        for(int i=1; i<ns; i++)
        {
            step(u,h,A[i],K[0],i,d,Z);
            effe(K[i],t+c[i]*h,Z);
            valf++;
        }

        step(u,h,b,K[0],ns,d,u);
        t+=h;
        //Se voglio stampare tutti gli stadi
        stampa(d,t,u,&prt);
        //stampa(d,t,u);
    }


    // inizio ciclo sul tempo
    for(unsigned long n=passi; n<=(N); n++)
    {
        copia(Old_U[(n-1)%passi],u,d);
        effe(Old_f[(n-1)%passi],t,u);
        valf++;

        Real coefu[passi];
        Real coeff[passi];

        for(int j=1; j<=passi; j++)
        {
            int k=(n-j)%passi;
            coeff[k]=sigma[j];
            coefu[k]=ro[j];
        }
        step(0,-1,coefu,Old_U[0],passi,d,u);
        step(u,h,coeff,Old_f[0],passi,d,u);
        t+=h;

        //Se voglio stampare tutti gli stadi
        stampa(d,t,u,&prt);
        //stampa(d,t,u);
    }

     //Se voglio stampare solo il dato iniziale e il dato finale
    //stampa(d,t,u,&prt);
        stampa(d,t,u);
    cout << "Numero valutazioni: " << valf << endl;
    return 0;
}
