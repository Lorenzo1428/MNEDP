clc
clear 
close all

dx = 0.05;
[u,X,Y,sol,eps,it] = fivePoints(dx);
err = max(max(abs(u - sol(X,Y))));

dx1 = 0.5*dx;
[u1,X1,Y1,sol1,eps1,it1] = fivePoints(dx1);
<<<<<<< HEAD
err1 = norm(u1 - sol1(X1,Y1),inf) + eps1;
S = sol(X1,Y1);
=======
err1 = max(max(abs(u1 - sol(X1,Y1))));

>>>>>>> 5f3505a (1.6)
log2(err/err1)

f1 = figure(Name="sol");
surf(X1,Y1,S);

f2 = figure(Name="approx");
surf(X1,Y1,u1);

function [u,X,Y,sol,eps,it] = fivePoints(dx)

    x = 0:dx:1;
    N = length(x);
    [X,Y] = ndgrid(x);
    
    f = @(x,y) (32*pi*pi*cos(4*pi*x).*sin(4*pi*y));
    phi = @(x) (sin(4*pi*x));
    sol = @(x,y) (cos(4*pi*x).*sin(4*pi*y));
    
    S = sol(X,Y);
    F = f(X,Y);
    u = zeros(N,N);
    
    u(1,:) = phi(x);
    u(:,1) = zeros(N,1);
    u(end,:) = phi(x);
    u(:,end) = zeros(N,1);
    
    it_max = 1e6;
    tol = 1e-8;
    it = 0;
    eps = 1;
    while eps > tol && it < it_max
        prev = u;
        for i = 2:N-1
            for j = 2:N-1
                u(i,j) = 0.25*(dx*dx*F(i,j) + u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1));
            end
        end
        it = it + 1;
        eps = max(max(abs(u - prev)));
    end

end