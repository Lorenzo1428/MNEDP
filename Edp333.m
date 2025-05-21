clc
clear
close all
E = EdpTrasportoD2;

dx = 0.5;
T = 0.5; 
a = [1 1];
[U,tf,dt,X,Y] = E.laxF2(dx,a,T);
sol = sin(X-tf*a(1)).*cos(Y-tf*a(2));
err1 = norm(U(1:end-1,1:end-1) - sol,inf);

dx1 = 0.5*dx;
[U,tf,dt,X,Y] = E.laxF2(dx1,a,T);
sol = sin(X-tf*a(1)).*cos(Y-tf*a(2));
err2 = norm(U(1:end-1,1:end-1) - sol,2);

p = log2(err1/err2);

f1 = figure(Name="Sol");
    surf(X,Y,sol);
f2 = figure;
    surf(X,Y,U(1:end-1,1:end-1));
    zlim([-1,1])

