clc
clear
close all

k = [1 2 3 4 5 6]
k(1:2:end)
k(2:2:end)


%% cond miste
dx = 0.0001;
x = 0:dx:1;
N = length(x);
f = @(x) 4*pi*pi*sin(2*pi*x);
sol = @(x) sin(2*pi*x);
epsi = 1;
beta = 2*pi;
F = dx*f(x(2:end)');
F(end) = 0.5*F(end) + beta;
A = (epsi/dx)*(2*eye(N-1) - diag(ones(1,N-2),-1) - diag(ones(1,N-2),1));
A(end,end) = 0.5*A(end,end);
U = [0;A\F];
%plot(x,U,x,sol(x))
err2 = sqrt(dx*sum(abs(U - sol(x'))));
errinf = max(abs(U - sol(x')));

%% neumann beta = 0, gamma non cost

dx = 0.1;
x = 0:dx:1;
xm = 0.5*(x(2:end) + x(1:end-1));
N = length(x);
f = @(x) sin(x);
epsi = 1;
gamma = @(x) 4 - x;
u0 = 2;
u1 = 1;
F = (1/3)*dx*(f(x(2:end-1))' + f(xm(1:end-1)') + f(xm(2:end)'));
F = [(1/6)*(f(x(1)) + 2*f(xm(1))) ; F; (1/6)*(f(x(end)) + 2*f(xm(end))) ];
F(1) = F(1) - u0;
F(end) = F(end) - u1;
A1 = (epsi/dx)*(2*eye(N) - diag(ones(1,N-1),-1) - diag(ones(1,N-1),1));
A1(1,1) = 0.5*A1(1,1);
A1(end,end) = 0.5*A1(end,end);

A2 = (dx/3)*(diag([0.5*gamma(x(1)) ,gamma(x(2:end-1)), 0.5*gamma(x(end))] + [0.5*gamma(xm(1)),0.5*gamma(xm(1:end-1)) + 0.5*gamma(xm(2:end)) ,0.5*gamma(xm(end))],0)  + diag(0.5*gamma(xm).*ones(1,N-1),-1) + diag(0.5*gamma(xm).*ones(1,N-1),1));
A = A1+A2;
U = A\F;



