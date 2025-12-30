% stimo errore finale e ordine di convergenza

clear all, close all
format long

p=4; %ordine del metodo AB4

%Carico i dati
x1=load('AB_1000');
x2=load('AB_2000');
x3=load('AB_4000');
x4=load('AB_8000');


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
coeff=2^p/(2^p-1);

%Calcolo gli errori finali
errfin1=norm(e2-e1)*coeff;
errfin2=norm(e3-e2)*coeff;
errfin3=norm(e4-e3)*coeff;

%Stimo l'ordine del metodo
p1=(log(errfin1)-log(errfin2))/log(2);
p2=(log(errfin2)-log(errfin3))/log(2);

