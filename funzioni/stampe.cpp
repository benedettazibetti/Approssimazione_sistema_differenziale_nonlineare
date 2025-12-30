#include <iostream>
#include <fstream>
#include "stampe.h"
using namespace std;
typedef double Real;
void stampamat(Real *a,int n,int m,ofstream *prt)
{
    if (prt==0)
    {
        for (int i=0; i<n; i++)
        {
            for (int j=0; j<m; j++)
            {
                cout << *(a+j+i*m) <<"\t";
            }
            cout << "\n";
        }
        cout<<" "<<endl;
    }
    else
    {
        for (int i=0; i<n; i++)
        {
            for (int j=0; j<m; j++)
            {
                *prt << *(a+j+i*m) <<"\t";
            }
            *prt << "\n";
        }
        *prt<< " "<<endl;
    }
}

void stampa(int d,Real t,Real *u,ofstream *prt)
{
    if (prt==0)
    {
        cout << t<< "\t";
        for(int i=0; i<d; i++)
        {
            cout << u[i]<<"\t";
        }
        cout <<endl;
    }
    else
    {
        *prt << t<< "\t";
        for(int i=0; i<d; i++)
        {
            *prt << u[i]<<"\t";
        }
        *prt <<endl;
    }
}

void copia (Real *X, Real *u, int d)
{
    for(int i=0; i<d; i++)
    {
        X[i]=u[i];
    }
}

