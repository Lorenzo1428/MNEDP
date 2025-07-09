clc
clear
close all

epsi = 0.00001;
sigma = 10;
mu = sqrt(sigma/epsi);
sol = @(x) (exp(mu*x) - exp(-mu*x))/(exp(mu) - exp(-mu));

dx = 0.01;
Pe = dx*dx*sigma/(epsi*6);
x = 0:dx:1;
N = length(x);
f = -dx*(sigma*x');
U = gal1(f,N,dx,epsi,sigma);
U = U + x';
f1 = figure;
plot(x,U);
title("Soluzione Pe = " + Pe);
legend("Soluzione Approssimata");
xlabel("x");
ylabel("U(x)");
xlim([-0,1.1]);
ylim([-0.1,1.1]);

U = gal1ML(f,N,dx,epsi,sigma);
U = U + x';
err1 = sqrt(dx*sum((U - sol(x')).^2));
f2 = figure;
plot(x,U);
title("Soluzione Mass-Lumping e Pe = " + Pe);
legend("Soluzione Approssimata");
xlabel("x");
ylabel("U(x)");
xlim([-0,1.1]);
ylim([-0.1,1.1]);

function U = gal1(f,N,dx,epsi,sigma)
    f = f(2:end-1);
    A1 = (epsi/dx)*(2*eye(N-2) - diag(ones(1,N-3),-1) - diag(ones(1,N-3),1));
    A2 = sigma*dx*((2/3)*eye(N-2) + (1/6)*diag(ones(1,N-3),-1) + (1/6)*diag(ones(1,N-3),1));
    A = A1 + A2;
    U = A\f;
    U = [0;U;0];
end

function U = gal1ML(f,N,dx,epsi,sigma)
    f = f(2:end-1);
    A1 = (epsi/dx)*(2*eye(N-2) - diag(ones(1,N-3),-1) - diag(ones(1,N-3),1));
    A2L = sigma*dx*eye(N-2);
    A = A1 + A2L;
    U = A\f;
    U = [0;U;0];
end