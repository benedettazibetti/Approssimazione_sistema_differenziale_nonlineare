%obbiettivo: disegnare le traiettorie dei 100 punti sulla sfera
clear all
close all

%Carico il file
F=load('G3_1000');

figure(1)
hold on
grid on

%Disegno la sfera di raggio 1 centrata nell'origine
[X,Y,Z]=sphere;
axis equal
surf(X,Y,Z,'FaceAlpha',0.2,'EdgeColor','none');

hold on
%Disegno gli 8 punti di equilibrio
plot3(1,0,0,"k*");
plot3(-1,0,0,"k*");

plot3(0,1,0,"k*");
plot3(0,-1,0,"k*");

plot3(0,0,1,"k*");
plot3(0,0,-1,"k*");

plot3(1/20,sqrt(399)/20,0,"k*");
plot3(1/20,-sqrt(399)/20,0,"k*");


%Disegno le traiettorie
k=0;
N=1000;
for i=1:100
    for j=1:N+1
        t(j)=F(k+j,1);
        x(j)=F(k+j,2);
        y(j)=F(k+j,3);
        z(j)=F(k+j,4);
    end
plot3(x,y,z);
k=k+N+1;
end