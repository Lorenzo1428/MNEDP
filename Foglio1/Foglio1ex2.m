clc
clear
close all
E = EdpTrasporto2D;

dx = 0.08;
T = 1; 
a = [1 1];
[U,dt,X,Y] = E.laxF2(dx,a,T);
sol = sin(X-T*a(1)).*cos(Y-T*a(2));
err = abs(U(1:end-1,1:end-1) - sol);
err1 = max(max(err));
err11 = dx*dx*sum(err,"all")

dx1 = 0.5*dx;
[U,dt,X,Y] = E.laxF2(dx1,a,T);
sol = sin(X-T*a(1)).*cos(Y-T*a(2));
err = abs(U(1:end-1,1:end-1) - sol);
err12 = dx1*dx1*sum(err,"all")
err2 = max(max(err));

p = log2(err1/err2);
p1 = log2(err11/err12);

f1 = figure(Name="Sol");
    surf(X,Y,sol);
f2 = figure;
    surf(X,Y,U(1:end-1,1:end-1));
    zlim([-1,1])
f3 = figure;
    mesh(err);

