clc
clear
close all

%% griglia unif

epsi = 1;
b = 1;
mu = b/epsi;
sol = @(x) (1 - exp(mu*x))./(1 - exp(mu));

dx = 0.01;
x = 0:dx:1;
N = length(x);
f = -b*dx*ones(N,1);
U = galerkin1_unif(f,dx,N,epsi,b);
U = U + x'; 
err = dx*sum((U - sol(x')).^2);
plot(x,U,x,sol(x'));

%% ordine 

dx = 0.5*dx;
x = 0:dx:1;
N = length(x);
f = -b*dx*ones(N,1);
U = galerkin1_unif(f,dx,N,epsi,b);
U = U + x'; 
err1 = dx*sum((U - sol(x')).^2);

p = log2(err/err1);

%% correzione upwind



function U = galerkin1_unif(f,dx,N,epsi,b)
    A1 = (epsi/dx)*(2*eye(N) - diag(ones(1,N-1),-1) - diag(ones(1,N-1),1));
    A2 = 0.5*b*(-diag(ones(1,N-1),-1) + diag(ones(1,N-1),1));
    A = A1 + A2;
    U = A\f;
end