% stimo errore finale e ordine di convergenza

clear all, close all
format long

p=4; %ordine del metodo RK4

%Carico i dati presi da CodeBlocks
x1=load('rk4_250');
x2=load('rk4_500');
x3=load('rk4_1000');
x4=load('rk4_2000');


%Calcolo gli errori
err1=x1(1,2:4)-x1(end,2:4);
e1=norm(err1);

err2=x2(1,2:4)-x2(end,2:4);
e2=norm(err2);

err3=x3(1,2:4)-x3(end,2:4);
e3=norm(err3);

err4=x4(1,2:4)-x4(end,2:4);
e4=norm(err4);

%Calcolo coefficiente di correzione (per compensazione errore sistematico)
coef=2^p/(2^p-1); 

%Calcolo gli errori finali
errfin1=norm(e2-e1)*coef;
errfin2=norm(e3-e2)*coef;
errfin3=norm(e4-e3)*coef;

%Stimo l'ordine del metodo 
p1=(log(errfin1)-log(errfin2))/log(2);
p2=(log(errfin2)-log(errfin3))/log(2);



