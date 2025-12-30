// PROGETTO ESAME N3
// METODO: RUNGE-KUTTA IMPLICITI
// Utilizzo GAUSS2   (ordine) p = 4     (numero di stadi) ns = 2


using namespace std;
#include <iostream>
#include <math.h>
#include <fstream>
#include "..\funzioni\mat_vet.h"
#include "..\funzioni\problemi.h"
#include "..\funzioni\newton.h"
#include "..\funzioni\stampe.h"
#include "..\funzioni\Runge_Kutta.h"

using namespace std;
typedef double Real;

const Real pi=4.0*atan(1.0);

void effeNewton(Real *F,Real *X);
void J_effeNewton(Real *J,Real *X);

// definisco puntatori al tipo di funzioni che enunciano il problema
void(*effe)(Real*,Real,Real*);
void(*dati)(Real *,Real *,Real *);
void(*butcher)(int, Real *,Real *,Real *,Real *, int *); //Carica i coefficienti della matrice di Butcher
void(*Jf)(Real*,Real,Real*);

// variabili comuni
Real t,T,h;

// dati del problema, del metodo e dimensioni
const int d=3;
Real u[d];
unsigned long valf;

const int ns=2;  //Numero degli stadi (ns=2 perchè sto considerando GAUSS2)
Real b[ns];
Real c[ns];
Real A[ns][ns];// matrice
Real K[ns][d];
//Real KK[ns][d];
Real Z[d]; // vettore

int main ()
{
    effe=f_progetto3;
    dati=dati_progetto3;
    Jf=Jf_progetto3;
    butcher=GAUSS2;
    butcher(ns,A[0],b,c,0,0);

    // predispongo il file stampa
    char n_file[21]= {0};
    cout << "dammi nome file di stampa(max 20 caratteri)";
    cin >> n_file ;
    ofstream prt(n_file);
    prt.precision(14);
    ofstream prt2("iterazioni");


    //parametri e variabili Newton
    Real toll=1e-14;
    int nitmax=20;
    valf=0;

    //input
    dati(&t,&T,u);
    unsigned long N;
    cout << "inserisci numero di passi N  = ";
    cin >> N;
    h=(T-t)/Real(N);


    for(int i=0;i<ns;i++){
        effe(K[i],t+c[i]*h,u);
    }

    // inizia ciclo sul tempo
    stampa(d,t,u);
    stampa(d,t,u,&prt);
    for(unsigned long n=1; n<=N; n++)
    {
        // iterazione Newton
        int nit =0;
        newton_sist(effeNewton,J_effeNewton, d*ns, K[0], &nit,toll,nitmax);
        step(u,h,b,K[0],ns,d,u);
        prt2<<t<<"\t"<<nit<<endl;
        // incremento il tempo
        t+=h;
        stampa(d,t,u,&prt);

    }

    stampa(d,t,u);

    cout << "valutazioni  f = "<< valf <<endl;
    return 0;
}

void effeNewton(Real *F,Real *KK)  //le scrivo nel main perchè sono tipiche di ogni metodo
{
    Real Z[d]={0};
    //Costruzione dei K (devo usare TUTTI i K)
    for(int i=0; i<ns; i++)
    {
        step(u,h,A[i],KK,ns,d,Z);
        effe(F+d*i,t+c[i]*h,Z);
        valf++;
    }
    for(int i=0; i<ns*d; i++)
    {
        F[i]=KK[i]-F[i]; //è la F corsivo
    }

}
void J_effeNewton(Real *J,Real *KK)
{
    Real Jff[d][d];
    Real Z[d]={0};
    //Calcolo il ciclo sulle righe
    for (int i=0; i<ns; i++)
    {
        step(u,h,A[i],KK,ns,d,Z);
        Jf(Jff[0],t+c[i]*h,Z);
        valf+=d;
        //Calcolo il ciclo sulle colonne
        for (int j=0; j<ns; j++)
        {
            Real aij= -A[i][j]*h;
            for(int l=0; l<d; l++)
            {
                for(int m=0; m<d; m++)
                {
                    int indice=(i*d+l)*(ns*d)+j*d+m;
                    J[indice]=aij*Jff[l][m];
                }
            }
        }
    }

    //J=J+I
    int k=0;
    for (int i=0; i<ns*d; i++)
    {
        J[k]++;
        k+=ns*d+1;
    }
}
