clc
clear 
close all

%% Diffusione Reazione 
epsi = 0.01;
sigma = 1;
mu = sqrt(sigma/epsi);
sol = @(x) (exp(mu*x) - exp(-mu*x))/(exp(mu) - exp(-mu));

dx = 0.001;
x = 0:dx:1;
N = length(x);
f = -dx*(sigma*x');
U = gal1(f,N,dx,epsi,sigma);
U = U + x';
err = sqrt(dx*sum((U - sol(x')).^2));
plot(x,U,x,sol(x'));
title("Soluzione con h= " + dx + ", eps= "+epsi + ", sigma= " +sigma );
legend("Soluzione approssimata","Soluzione esatta");
xlabel("x");
ylabel("U(x)");
xlim([-0,1.1]);
ylim([-0.1,1.1]);

%% Ordine
M = 8;
dx1 = zeros(M,1);
p = zeros(M,1);
dx1(1) = 0.1;
err1 = zeros(M,1);
for i = 1:M
    x = 0:dx1(i):1;
    N = length(x);
    f = -dx1(i)*(sigma*x');
    U = gal1(f,N,dx1(i),epsi,sigma);
    U = U + x';
    err1(i) = sqrt(dx1(i)*sum((U - sol(x')).^2));
    if i>1
        p(i) = log2(err1(i-1)/err1(i));
    end
    if i < M
        dx1(i+1) = 0.5*dx1(i);
    end
end
table(dx1,err1,p)

dx1ML = zeros(M,1);
pML = zeros(M,1);
dx1ML(1) = 0.1;
err1ML = zeros(M,1);
for i = 1:M
    x = 0:dx1ML(i):1;
    N = length(x);
    f = -dx1ML(i)*(sigma*x');
    U = gal1ML(f,N,dx1ML(i),epsi,sigma);
    U = U + x';
    err1ML(i) = sqrt(dx1ML(i)*sum((U - sol(x')).^2));
    if i>1
        pML(i) = log2(err1ML(i-1)/err1ML(i));
    end
    if i < M
        dx1ML(i+1) = 0.5*dx1ML(i);
    end
end
table(dx1ML,err1ML,pML)

f2 = figure;
semilogy(1:M,err1,1:M,err1ML)
title("Andamento errore in L2");
legend("Galerkin","Mass-lumping");


%% galerkin dim 1
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