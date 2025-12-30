using namespace std;
typedef double Real;
#include <iostream>
#include <string.h>
#include "..\funzioni\multi_step_metodi.h"
void normalizza(Real *ro,Real *sigma,int passi);
//
void AB(int ordine,Real ro[],Real sigma[]){

    // AB sta per Adams_Bashforth, ESPLICITI
    // l'intero in input è l'ordine del metodo
    // il numero dei passi coincide con l'ordine
    if((ordine<1)||(ordine>6)){
        cout<<"metodo non previsto, processo abortito"<<endl;
        return;
    }
    switch(ordine){
    case 1:{
                ro[0]=1.0;
                ro[1]=-1.0;
                sigma[0]=0.0;
                sigma[1]=1.0;
                }
                break;
    case 2:{
                ro[0]=1.0;
                ro[1]=-1.0;
                ro[2]=0.0;
                sigma[0]=0.0;
                sigma[1]=3.0/2.0;
                sigma[2]=-0.5;
                }
                break;
    case 3:{
                ro[0]=12.0;
                ro[1]=-12.0;
                ro[2]=0.0;
                ro[3]=0.0;
                sigma[0]=0.0;
                sigma[1]=23.0;
                sigma[2]=-16.0;
                sigma[3]=5.0;
                }
                break;
    case 4:{
                ro[0]=24.0;
                ro[1]=-24.0;
                ro[2]=0.0;
                ro[3]=0.0;
                ro[4]=0.0;
                sigma[0]=0.0;
                sigma[1]=55.0;
                sigma[2]=-59.0;
                sigma[3]=37.0;
                sigma[4]=-9.0;
                }
                break;
    case 5:{
                ro[0]=720.0;
                ro[1]=-720.0;
                ro[2]=0.0;
                ro[3]=0.0;
                ro[4]=0.0;
                ro[5]=0.0;
                sigma[0]=0.0;
                sigma[1]=1901.0;
                sigma[2]=-2774.0;
                sigma[3]=2616.0;
                sigma[4]=-1274;
                sigma[5]=251.0;
                }
                break;
    case 6:{
                ro[0]=1440.0;
                ro[1]=-1440.0;
                ro[2]=0.0;
                ro[3]=0.0;
                ro[4]=0.0;
                ro[5]=0.0;
                ro[6]=0.0;
                sigma[0]=0.0;
                sigma[1]=4277.0;
                sigma[2]=-7923.0;
                sigma[3]=9982.0;
                sigma[4]=-7298;
                sigma[5]=2877.0;
                sigma[6]=-475.0;
                }
                break;
    }
   normalizza(ro,sigma,ordine);
}

void AM(int ordine,Real ro[],Real sigma[]){

    // AM sta per Adams_Moulton, IMPLICITI
    // l'intero in input è l'ordine del metodo
    // il numero dei passi = ordine -1;
    if((ordine<2)||(ordine>4)){
        cout<<"metodo non previsto, processo abortito"<<endl;
        return;
    }
    switch(ordine){
    case 2:{
                ro[0]=1.0;
                ro[1]=-1.0;
                sigma[0]=0.5;
                sigma[1]=0.5;
                }
                break;
            case 3:{
                ro[0]=12.0;
                ro[1]=-12.0;
                ro[2]=0.0;
                sigma[0]=5.0;
                sigma[1]=8.0;
                sigma[2]=-1.0;
                }
                break;
            case 4:{
                ro[0]=24.0;
                ro[1]=-24.0;
                ro[2]=0.0;
                ro[3]=0.0;
                sigma[0]=9.0;
                sigma[1]=19.0;
                sigma[2]=-5.0;
                sigma[3]=1.0;
                }
                break;
    }
    normalizza(ro,sigma,ordine-1);
  }
void BDF(int ordine,Real ro[],Real sigma[]){

    // BDF sta per Backward Differentiation formula IMPLICITI
    // l'intero in input è l'ordine del metodo
    // il numero dei passi coincide con l'ordine;
    if((ordine<2)||(ordine>6)){
        cout<<"metodo non previsto, processo abortito"<<endl;
        return;
    }
    switch(ordine){
    case 2:{
                ro[0]=1.5;
                ro[1]=-2.0;
                ro[2]=0.5;
                sigma[0]=1.0;
                sigma[1]=0.0;
                sigma[2]=0.0;
                }
                break;
            case 3:{
                ro[0]=11.0/6.0;
                ro[1]=-3.0;
                ro[2]=1.5;
                ro[3]=-1.0/3.0;
                sigma[0]=1.0;
                sigma[1]=0.0;
                sigma[2]=0.0;
                sigma[3]=0.0;
                }
                break;
            case 4:{
                ro[0]=25.0/12.0;
                ro[1]=-4.0;
                ro[2]=3.0;
                ro[3]=-4.0/3.0;
                ro[4]=0.25;
                sigma[0]=1.0;
                sigma[1]=0.0;
                sigma[2]=0.0;
                sigma[3]=0.0;
                sigma[4]=0.0;
                }
                break;
  case 5:{
                ro[0]=137.0/60.0;
                ro[1]=-5.0;
                ro[2]=5.0;
                ro[3]=-10.0/3.0;
                ro[4]=1.25;
                ro[5]=-0.2;
                sigma[0]=1.0;
                sigma[1]=0.0;
                sigma[2]=0.0;
                sigma[3]=0.0;
                sigma[4]=0.0;
                sigma[5]=0.0;
                }
                break;
            case 6:   {
				ro[0]=49.0/20.0;
                ro[1]=-6.0;
                ro[2]=15.0/2.0;
                ro[3]=-20.0/3.0;
                ro[4]=15.0/4.0;
                ro[5]=-6.0/5.0;
                ro[6]=1.0/6.0;
                sigma[0]=1.0;
                sigma[1]=0.0;
                sigma[2]=0.0;
                sigma[3]=0.0;
                sigma[4]=0.0;
                sigma[5]=0.0;
                sigma[6]=0.0;
                }
                break;
                }
    normalizza(ro,sigma,ordine);
}

void LF(int ordine,Real ro[],Real sigma[]){

    // LF sta per Leap Frog ESPLICITO
    // l'intero in input è l'ordine del metodo
    // il numero dei passi coincide con l'ordine;
                ro[0]=1.0;
                ro[1]=0.0;
                ro[2]=-1.0;
                sigma[0]=0.0;
                sigma[1]=2.0;
                sigma[2]=0.0;
}
void normalizza(Real *ro,Real *sigma,int passi)
{
    if(ro[0]!=1.0){
        Real aux=ro[0];
            for(int i=0;i<=passi;i++){
                ro[i]=ro[i]/aux;
                sigma[i]=sigma[i]/aux;
            }
       ro[0]=1.0;
    }
}
