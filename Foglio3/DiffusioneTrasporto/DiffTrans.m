clc
clear
close all

%% Diffusione Trasporto

epsi = 0.1;
b = 1;
mu = b/epsi;
sol = @(x) (exp(mu*x)-1)./(exp(mu)-1);

dx = 0.001;
x = 0:dx:1;
N = length(x);
f = -b*dx*ones(N-2,1);
U = galerkin1_unif(f,dx,N,epsi,b);
U = U + x'; 
err = sqrt(dx*sum((U - sol(x')).^2));

%griglia non uniforme
nodi = x.*exp(1-x);
h = nodi(2:end)-nodi(1:end-1);
hmax = max(h);
f = -0.5*b*(h(1:end-1) + h(2:end))';
Uh = galerkin1_non_unif(f,h,N,epsi,b);
Uh = Uh + nodi'; 
errh = sqrt(hmax*sum((U - sol(nodi')).^2));

f1 = figure;
subplot(2,1,1)
plot(nodi,Uh,nodi,sol(nodi'));
legend("soluzione approssimata","sol esatta");
title("Soluzione Galerkin con griglia non uniforme")
xlim([0,1.1]);
ylim([-0.1,1.1]);

subplot(2,1,2)
plot(x,U,x,sol(x'));
legend("soluzione approssimata","sol esatta");
title("Soluzione Galerkin con griglia uniforme")
xlim([0,1.1]);
ylim([-0.1,1.1]);

%% ordine 

M = 8;
dx1 = zeros(M,1);
p = zeros(M,1);
dx1(1) = 0.1;
err1 = zeros(M,1);
for i = 1:M
    x = 0:dx1(i):1;
    N = length(x);
    f = -b*dx1(i)*ones(N-2,1);
    U = galerkin1_unif(f,dx1(i),N,epsi,b);
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

dx1Up = zeros(M,1);
pUp = zeros(M,1);
dx1Up(1) = 0.1;
epsil = epsi + 0.5*b*dx1Up(1);
err1Up = zeros(M,1);
for i = 1:M
    x = 0:dx1Up(i):1;
    N = length(x);
    f = -b*dx1Up(i)*ones(N-2,1);
    U = galerkin1_unif(f,dx1Up(i),N,epsil,b);
    U = U + x'; 
    err1Up(i) = sqrt(dx1(i)*sum((U - sol(x')).^2));
    if i>1
        pUp(i) = log2(err1Up(i-1)/err1Up(i));
    end
    if i < M
        dx1Up(i+1) = 0.5*dx1Up(i);
        epsil = epsi + 0.5*b*dx1Up(i);
    end
end
table(dx1Up,err1Up,pUp)

dx1N = zeros(M,1);
pN = zeros(M,1);
dx1N(1) = 0.1;
err1N = zeros(M,1);
hmax1 = zeros(M,1);
for i = 1:M
    x = 0:dx1N(i):1;
    N = length(x);
    nodi1 = x.*exp(1-x);
    h1 = nodi1(2:end)-nodi1(1:end-1);
    hmax1(i) = max(h1);
    f = -0.5*b*(h1(1:end-1) + h1(2:end))';
    Uh = galerkin1_non_unif(f,h1,N,epsi,b);
    Uh = Uh + nodi1'; 
    err1N(i) = sqrt(hmax1(i)*sum((Uh - sol(nodi1')).^2));
    if i>1
        pN(i) = log2(err1N(i-1)/err1N(i));
    end
    if i < M
        dx1N(i+1) = 0.5*dx1N(i);
    end
end
table(hmax1,err1N,pN)

f2 = figure;
semilogy(1:M,err1,1:M,err1Up,1:M,err1N)
title("Andamento errore in L2");
legend("Galerkin","Stabilizazzione Upwind","Galerkin non uniforme");

%% Galerkin

function U = galerkin1_unif(f,dx,N,epsi,b)
    A1 = (epsi/dx)*(2*eye(N-2) - diag(ones(1,N-3),-1) - diag(ones(1,N-3),1));
    A2 = 0.5*b*(-diag(ones(1,N-3),-1) + diag(ones(1,N-3),1));
    A = A1 + A2;
    U = A\f;
    U = [0;U;0];
end

function U = galerkin1_non_unif(f,h,N,epsi,b)
    B = [1 -1;-1 1];
    B1 = zeros(N);
    for i = 1:N-1
        B1(i:i+1,i:i+1) = B1(i:i+1,i:i+1) + B/h(i);
    end
    A1 = epsi*B1(2:end-1,2:end-1);
    A2 = 0.5*b*(-diag(ones(1,N-3),-1) + diag(ones(1,N-3),1));
    A = A1 + A2;
    U = A\f;
    U = [0;U;0];
end