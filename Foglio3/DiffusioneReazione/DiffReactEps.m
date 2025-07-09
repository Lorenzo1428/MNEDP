clc
clear 
close all

epsi = 0.01;
sigma = 1;
mu = sqrt(sigma/epsi);

sol = @(x) (exp(mu*x) - exp(-mu*x))/(exp(mu) - exp(-mu));

dx = 0.001;
x = 0:dx:1;
N = length(x);
f = -dx*(sigma*x');
epsi_ = 10;
U_ = zeros(N,6);
for i = 1:6
    epsi_ = 0.1*epsi_;
    U_(:,i) = gal1(f,N,dx,epsi_,sigma);
    U_(:,i) = U_(:,i) + x';
end
plot(x,U_)
title("Soluzione approssimata al variare di eps")
legend("eps = 1","eps = 0.1", "eps = 0.01", "eps = 1e-3", "eps = 1e-4", "eps = 1e-5");
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