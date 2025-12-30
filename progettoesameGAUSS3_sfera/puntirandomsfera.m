%estraggo 100 punti casuali sulla superficie sferica
clear all, close all
format long;

rng(0,'twister');             %Inizializzo il generatore di numeri random
rvals = 2*rand(100,1)-1;      % vettore random 100x1 
elevation = asin(rvals);      %Seno alla meno 1 del vettore
azimuth = 2*pi*rand(100,1);   %Angolo azimutale
%radii=1*rand(100,1).^(1/3);  %Valore del raggio per ogni punto della sfera

%Converto in coordinate cartesiane e disegno i punti sulla sfera
[x,y,z] = sph2cart(azimuth,elevation,1);
figure(1)
plot3(x,y,z,'.')
grid on
axis equal
%Disegno la sfera di raggio 1 centrata nell'origine
hold on
[X,Y,Z]=sphere;
axis equal
surf(X,Y,Z,'FaceAlpha',0.2,'EdgeColor','none');
hold off

%Salvo i punti sul file che mi serviranno in CodeBlocks
fileid= fopen('puntirandomsfera.txt','wt'); %'wt' per aprire file con TextEditor

for i=1:100
    X(i)=x(i);
    Y(i)=y(i);
    Z(i)=z(i);
    fprintf(fileid,'%16.16g\t %16.16g\t %16.16g\n', [X(i),Y(i),Z(i)]); %The more compact of %e or %f, with no trailing zeros 
end

fclose(fileid);