clc
clear
close all

dx = 0.1;
epsi = 0.1;
b = 10;
mu = b/epsi;
x = 0:dx:1;
N = length(x);
f = -b*dx*ones(N-2,1);
sol = @(x) (exp(mu*x)-1)./(exp(mu)-1);
Pe=0.5*dx*mu;
U = galerkin1_unif(f,dx,N,epsi,b);
U = U + x'; 
err = sqrt(dx*sum((U - sol(x')).^2));
f1 = figure;
plot(x,U,x,sol(x'));
title("Soluzione galerkin con Pe= " + Pe);
legend("soluzione approssimata","sol esatta");
xlabel("x");
ylabel("U(x)");
xlim([-0,1.1]);
ylim([-0.1,1.1]);


epsil = epsi + 0.5*b*dx;
mul = b/epsil;
Pe1=0.5*dx*mul;
sol1 = @(x) (exp(mul*x)-1)./(exp(mul)-1);
Up = galerkin1_unif(f,dx,N,epsil,b);
Up = Up + x'; 
errUp = sqrt(dx*sum((U - sol1(x')).^2));
f2 = figure;
plot(x,Up,x,sol1(x'));
title("Correzione Upwind con Pe^{up}= "+ Pe1);
legend("soluzione approssimata","sol esatta");
xlabel("x");
ylabel("U(x)");
xlim([-0,1.1]);
ylim([-0.1,1.1]);

function U = galerkin1_unif(f,dx,N,epsi,b)
    A1 = (epsi/dx)*(2*eye(N-2) - diag(ones(1,N-3),-1) - diag(ones(1,N-3),1));
    A2 = 0.5*b*(-diag(ones(1,N-3),-1) + diag(ones(1,N-3),1));
    A = A1 + A2;
    U = A\f;
    U = [0;U;0];
end