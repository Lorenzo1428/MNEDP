clc
clear 
close all

id = 3;

switch id
    case 1
        a = 0;
        b = 1;
        f = @(x,y) zeros(length(x),length(y));
        sol = @(x,y) exp(pi*x).*(cos(pi*y));
        phi = @(x,y) (exp(pi*x)).*(y == 0 & x ~= 0) + (-exp(pi*x)).*(y == 1 & x ~= 1) ...
            + (cos(pi*y)).*(x == 0 & y ~= 1) + (exp(pi)*cos(pi*y)).*(x == 1 & y ~= 0);

    case 2
        a = 0;
        b = 1;
        f = @(x,y) (32*pi*pi*cos(4*pi*x).*sin(4*pi*y));
        phi = @(x,y) (sin(4*pi*y));
        sol = @(x,y) (cos(4*pi*x).*sin(4*pi*y));
    case 3
        a = -1;
        b = 1;
        f = @(x,y) zeros(length(x),length(y));
        sol = @(x,y) 2*(1 + y)./((3 + x).^2 + (1 + y).^2);
        phi = @(x,y) zeros(length(x)).*(y == -1 & x ~= -1) + 4./((3 + x).^2 + 4).*(y == 1 & x ~= 1) ...
            + 2*(1 + y)./(4 + (1 + y).^2).*(x == -1 & y ~= 1) + 2*(1 + y)./(16 + (1 + y).^2).*(x == 1 & y ~= -1);
end

dx = 0.2;
N = floor((b-a)/dx);
dx = (b-a)/N;
x = a:dx:b;
[X,Y] = ndgrid(x);
N = length(x);
[u,it] = ninePoints(dx,X,Y,N,f,phi);
err = max(max(abs(u - sol(X,Y))));

dx1 = 0.5*dx;
N = floor((b-a)/dx);
dx = (b-a)/N;
x = a:dx1:b;
[X1,Y1] = ndgrid(x);
N = length(x);
[u1,it1] = ninePoints(dx1,X1,Y1,N,f,phi);

S = sol(X1,Y1);
err1 = max(max(abs(u1 - sol(X1,Y1))));
log2(err/err1)

f1 = figure(Name="sol");
surf(X1,Y1,S);
title("Soluzione esatta");
f2 = figure(Name="approx");
surf(X1,Y1,u1);
title("Soluzione approssimata");