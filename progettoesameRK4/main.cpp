// PROGETTO ESAME N3
// METODO: RUNGE-KUTTA ESPLICITO
// Utilizzo RK4 (ordine) p=4    (numero di stadi) ns=4

#include <iostream>
#include <math.h>
#include <fstream>
#include "..\funzioni\problemi.h"
#include "..\funzioni\mat_vet.h"
#include "..\funzioni\Runge_Kutta.h"
#include "..\funzioni\stampe.h"

using namespace std;
typedef double Real;
const Real pi=4.0*atan(1.0);

int valf=0;

// definisco puntatori al tipo di funzioni che enunciano il problema
void(*effe)(Real*, Real,Real*); //t, y, f(t)
void(*dati)(Real *,Real *,Real *);
void(*butcher)(int, Real *,Real *,Real *,Real *, int *); //Carica i coefficienti della matrice di Butcher

// variabili comuni
Real t,T,h;


// dati del problema e dimensioni
const int d=3; // sfera:problema tridimensionale
Real u[d]= {0};

const int ns=4; //Numero degli stadi (ns=4 perchè sto considerando RK4)
Real b[ns];
Real c[ns];
Real A[ns][ns];
Real K[ns][d];
Real Z[d];

int main ()
{
    effe=f_progetto3;
    dati=dati_progetto3;
    butcher=RK4;
    butcher(ns,A[0],b,c,0,0);

    // predispongo il file stampa
    char n_file[21]= {0};
    cout << "dammi nome file di stampa(max 20 caratteri)";
    cin >> n_file ;
    ofstream prt(n_file);
    prt.precision(14);


    // input
    dati(&t,&T,u);
    unsigned long N; //Numero intero privo di segno
    cout << "inserisci numero di passi N  = ";
    cin >> N;
    h=(T-t)/Real(N); //Real perchè non divido per interi ma divido per reali

    stampa(d,t,u,&prt);
    stampa(d,t,u);



    // inizio ciclo sul tempo

    for(unsigned long n=1; n<=N; n++)
    {

        for(int i=0; i<ns; i++)
        {
            step(u,h,A[i],K[0],i,d,Z);
            effe(K[i],t+c[i]*h,Z);
            valf=valf+1;
        }

        step(u,h,b,K[0],ns,d,u);
        t+=h;

        //Se voglio stampare tutti i passaggi intermedi
        stampa(d,t,u,&prt);
        //stampa(d,t,u);
    }
        //Se voglio stampare solo i dati finali
        //stampa(d,t,u,&prt);
        stampa(d,t,u);

        cout << "numero valutazioni f  = " << valf << endl;


       return 0;
}
