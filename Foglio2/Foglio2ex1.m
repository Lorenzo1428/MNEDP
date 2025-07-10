clc
clear 
close all

disp("1. punto a. 2. punto b 3. punto c")
prompt = "Scegli un problema: ";
id = input(prompt);

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

dx = 0.05;
N = floor((b-a)/dx);
dx = (b-a)/N;
x = a:dx:b;
N = length(x);
[X,Y] = ndgrid(x);
S = sol(X,Y);
[u5,it5] = fivePoints(dx,X,Y,N,f,phi);
[u9,it9] = ninePoints(dx,X,Y,N,f,phi);
err5 = max(max(abs(u5 - S)));
err9 = max(max(abs(u9 - S)));

f1 = figure(Name="Grafici soluzioni approssimate");
subplot(1,2,1)
surf(X,Y,u5);
title("Lapaciano 5 punti");
subplot(1,2,2)
surf(X,Y,u5);
title("Lapaciano 9 punti");

f2 = figure(Name="Soluzione esatta e errori");
surf(X,Y,S);
title("Soluzione esatta");

f3 = figure(Name="Errore in L2");
subplot(1,2,1)
mesh(abs(u5 - S));
title("Errore laplaciano 5 punti");
subplot(1,2,2)
mesh(abs(u9 - S))
title("Errore laplaciano 9 punti");