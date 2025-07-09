clc
clear
close all

b = 1;
dx = 0.0001;
x = 0:dx:1;
N = length(x);
f = -b*dx*ones(N-2,1);

U = zeros(N,6);
epsi = 10;
for i = 1:6
    epsi = 0.1*epsi;
    mu = b/epsi;
    Pl=0.5*dx*mu;
    U(:,i) = galerkin1_unif(f,dx,N,epsi,b);
    U(:,i) = U(:,i) + x'; 
end

f1 = figure;
plot(x,U);
title("Soluzione Galerkin");
legend("eps = 1","eps = 0.1","eps = 0.01","eps = 0.001","eps = 10^{-4}","eps = 10^{-5}");
xlim([0,1.1]);
ylim([-0.1,1.1]);

epsi = 10;
epsil = epsi + 0.5*b*dx;
for i = 1:6
    epsi = 0.1*epsi;
    epsil = epsi + 0.5*b*dx;
    mu = b/epsi;
    Pl=0.5*dx*mu;
    U(:,i) = galerkin1_unif(f,dx,N,epsil,b);
    U(:,i) = U(:,i) + x'; 
end

f2 = figure;
plot(x,U);
title("Soluzione Upwind");
legend("eps = 1","eps = 0.1","eps = 0.01","eps = 0.001","eps = 10^{-4}","eps = 10^{-5}");
xlim([0,1.1]);
ylim([-0.1,1.1]);



function U = galerkin1_unif(f,dx,N,epsi,b)
    A1 = (epsi/dx)*(2*eye(N-2) - diag(ones(1,N-3),-1) - diag(ones(1,N-3),1));
    A2 = 0.5*b*(-diag(ones(1,N-3),-1) + diag(ones(1,N-3),1));
    A = A1 + A2;
    U = A\f;
    U = [0;U;0];
end