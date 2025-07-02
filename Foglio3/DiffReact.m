clc
clear 
close all

%% dati iniziali 
a = 0;
b = 1;
ua = 0;
ub = 1;
epsi = 0.01;
sigma = 1;
mu = sqrt(sigma/epsi);

r = @(x) ua*(b - x)/(b-a) + ub*(x - a)/(b-a);
sol = @(x) (exp(mu*x) - exp(-mu*x))/(exp(mu) - exp(-mu));

%% eps = ...
dx = 0.001;
x = a:dx:b;
N = length(x);
f = -dx*(sigma*r(x'));
epsi_ = 10;
U_ = zeros(N,6);
for i = 1:6
    epsi_ = 0.1*epsi_;
    U_(:,i) = gal1(f,N,dx,epsi_,sigma);
    U_(:,i) = U_(:,i) + r(x');
end
plot(x,U_)
legend("eps = 1","eps = 0.1", "eps = 0.01", "eps = 1e-3", "eps = 1e-4", "eps = 1e-5");

%% ordine in L2
% dx
dx = 0.001;
x = a:dx:b;
N = length(x);
f = -dx*(sigma*r(x'));
U = gal1(f,N,dx,epsi,sigma);
U = U + r(x');
err = dx*sum((U - sol(x')).^2);
%err = norm(U - sol(x'), inf);
%plot(x,U,x,sol(x'))

% 0.5dx
dx = 0.5*dx;
x = a:dx:b;
N = length(x);
f = -dx*(sigma*r(x'));
U = gal1(f,N,dx,epsi,sigma);
U = U + r(x');
err1 = dx*sum((U - sol(x')).^2);
%err1 = norm(U - sol(x'), inf);

% ordine
p = log2(err/err1);


%% galerkin dim 1
function U = gal1(f,N,dx,epsi,sigma)
    A1 = (epsi/dx)*(2*eye(N) - diag(ones(1,N-1),-1) - diag(ones(1,N-1),1));
    A2 = sigma*dx*((2/3)*eye(N) + (1/6)*diag(ones(1,N-1),-1) + (1/6)*diag(ones(1,N-1),1));
    A = A1 + A2;
    U = A\f;
end