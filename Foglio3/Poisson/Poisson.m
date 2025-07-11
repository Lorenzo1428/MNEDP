clc
clear
close all

F = @(x) (0.25*pi^2)*sin(0.5*pi*x);
sol = @(x) sin(0.5*pi*x);

%% uniforme lineare

dx = 0.001;
x = 0:dx:1;
N = length(x);
f = F(x')*dx;
A = (1/dx)*(2*eye(N-2) - diag(ones(1,N-3),-1) - diag(ones(1,N-3),1));
U = [0; A\f(2:end-1) ; 0] + x';
err2_dx = sqrt(dx*sum((U - sol(x')).^2));

%% non uniforme lineare random

x1 = [0; sort(rand(N-2,1)); 1];
h = x1(2:end) - x1(1:end-1);
f1 = 0.5*F(x1(2:end-1)).*(h(1:end-1) + h(2:end));
B = [1 -1;-1 1];
B1 = zeros(N);
for i = 1:N-1
    B1(i:i+1,i:i+1) = B1(i:i+1,i:i+1) + B/h(i);
end
A1 = B1(2:end-1,2:end-1); 
U1 = A1\f1;
U1 = [0; U1; 0] + x1;
err2_h = sqrt(max(h)*sum((U1 - sol(x1)).^2));

%% elementi quadratici

f = zeros(N-2,1);
f(1:2:end) = (4/3)*F(x(2:2:end-1)')*dx;
f(2:2:end) = (2/3)*F(x(3:2:end-1)')*dx;
C0 = 1/3*[7 -8 1; -8 16 -8; 1 -8 7];
C = zeros(N);
for i = 1:2:N-2
    C(i:i+2,i:i+2) = C(i:i+2,i:i+2) + 0.5*C0/dx;
end
Aq = C(2:end-1,2:end-1);
Uq = [0; Aq\f; 0] + x';
err2_dxq = sqrt(dx*sum((Uq - sol(x')).^2));

f1 = figure;
subplot(3,1,1)
plot(x,U,x,sol(x'));
title("Elementi finiti lineari e griglia uniforme")
legend("Soluzione approssimata","Soluzione esatta");

subplot(3,1,2)
plot(x,U1,x,sol(x1'));
title("Elementi finiti lineari e griglia non uniforme")
legend("Soluzione approssimata","Soluzione esatta");

subplot(3,1,3)
plot(x,Uq,x,sol(x'));
title("Elementi finiti quadratici e griglia uniforme")
legend("Soluzione approssimata","Soluzione esatta");


%% ordine uniforme lineare
M = 7;

dx1 = zeros(M,1);
p = zeros(M,1);
dx1(1) = 0.1;
err1 = zeros(M,1);
for i = 1:M
    x = 0:dx1(i):1;
    N = length(x);
    f = F(x')*dx1(i);
    A = (1/dx1(i))*(2*eye(N-2) - diag(ones(1,N-3),-1) - diag(ones(1,N-3),1));
    U = [0; A\f(2:end-1) ; 0] + x';
    err1(i) = sqrt(dx1(i)*sum((U - sol(x')).^2));
    if i>1
        p(i) = log2(err1(i-1)/err1(i));
    end
    if i < M
        dx1(i+1) = 0.5*dx1(i);
    end
end
table(dx1,err1,p)

%% ordine quadratico 

dx1q = zeros(M,1);
pq = zeros(M,1);
dx1q(1) = 0.1;
err1q = zeros(M,1);
for i = 1:M
    x = 0:dx1q(i):1;
    N = length(x);
    f = zeros(N-2,1);
    f(1:2:end) = (4/3)*F(x(2:2:end-1)')*dx1q(i);
    f(2:2:end) = (2/3)*F(x(3:2:end-1)')*dx1q(i);
    C0 = 1/3*[7 -8 1; -8 16 -8; 1 -8 7];
    C = zeros(N);
    for j = 1:2:N-2
        C(j:j+2,j:j+2) = C(j:j+2,j:j+2) + 0.5*C0/dx1q(i);
    end
    Aq = C(2:end-1,2:end-1);
    U = [0; Aq\f; 0] + x';
    err1q(i) = sqrt(dx1q(i)*sum((U - sol(x')).^2));
    if i>1
        pq(i) = log2(err1q(i-1)/err1q(i));
    end
    if i < M
        dx1q(i+1) = 0.5*dx1q(i);
    end
end
table(dx1q,err1q,pq)

%% ordine non uniforme

dx1n = zeros(M,1);
pn = zeros(M,1);
dx1n(1) = 0.1;
err1n = zeros(M,1);
for i = 1:M
    xn = 0:dx1n(i):1;
    N = length(xn);
    x1n = [0; sort(rand(N-2,1)); 1];
    h = x1n(2:end) - x1n(1:end-1);
    f1 = 0.5*F(x1n(2:end-1)).*(h(1:end-1) + h(2:end));
    B = [1 -1;-1 1];
    B1 = zeros(N);
    for j = 1:N-1
        B1(j:j+1,j:j+1) = B1(j:j+1,j:j+1) + B/h(j);
    end
    A1 = B1(2:end-1,2:end-1); 
    U1 = A1\f1;
    U1 = [0; U1; 0] + x1n;
    err1n(i) = sqrt(max(h)*sum((U1 - sol(x1n)).^2));
    if i>1
        pn(i) = log2(err1n(i-1)/err1n(i));
    end
    if i < M
        dx1n(i+1) = 0.5*dx1n(i);
    end
end
table(dx1n,err1n,pn)


f2 = figure;
semilogy(1:M,err1,1:M,err1q,1:M,err1n)
title("Andamento errore in L2");
legend("Elementi lineari","Elementi quadratici","Lineari non uniformi");