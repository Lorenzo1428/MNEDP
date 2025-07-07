clc
clear
close all

F = @(x) (0.25*pi^2)*sin(0.5*pi*x);
sol = @(x) sin(0.5*pi*x);

%% uniforme lineare

dx = 0.002;
x = 0:dx:1;
N = length(x);
f = F(x')*dx;
A = (1/dx)*(2*eye(N-2) - diag(ones(1,N-3),-1) - diag(ones(1,N-3),1));
U = [0; A\f(2:end-1) ; 0] + x';
err2_dx = sqrt(dx*sum((U - sol(x')).^2));

%% non uniforme lineare random

x1 = [0; sort(rand(N-2,1)); 1];%x'.*exp(1-x');
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

%% ordine uniforme lineare

dx = 0.5*dx;
x = 0:dx:1;
N = length(x);
f = F(x')*dx;
A = (1/dx)*(2*eye(N-2) - diag(ones(1,N-3),-1) - diag(ones(1,N-3),1));
U = [0; A\f(2:end-1) ; 0] + x';
err2_dxhalf = sqrt(dx*sum((U - sol(x')).^2));
err1_dxhalf = dx*sum(abs(U - sol(x')));
p2 = log2(err2_dx/err2_dxhalf);

%% elementi quadratici

dx = 0.005;
x = 0:dx:1;
N = length(x);
f = zeros(N-2,1);
f(1:2:end) = (4/3)*F(x(2:2:end-1)')*dx;
f(2:2:end) = (2/3)*F(x(3:2:end-1)')*dx;
C0 = 1/3*[7 -8 1; -8 16 -8; 1 -8 7];
C = zeros(N);
for i = 1:2:N-2
    C(i:i+2,i:i+2) = C(i:i+2,i:i+2) + 0.5*C0/dx;
end
Aq = C(2:end-1,2:end-1);
U = [0; Aq\f; 0] + x';
err2_dxq = sqrt(dx*sum((U - sol(x')).^2));
plot(x,U,x,sol(x'));
legend("approx","sol");

%% ordine quadratico 

dx = 0.5*dx;
x = 0:dx:1;
N = length(x);
f = zeros(N-2,1);
f(1:2:end) = (4/3)*F(x(2:2:end-1)')*dx;
f(2:2:end) = (2/3)*F(x(3:2:end-1)')*dx;
C0 = 1/3*[7 -8 1; -8 16 -8; 1 -8 7];
C = zeros(N);
for i = 1:2:N-2
    C(i:i+2,i:i+2) = C(i:i+2,i:i+2) + 0.5*C0/dx;
end
Aq = C(2:end-1,2:end-1);
U = [0; Aq\f; 0] + x';
err2_dxhalfq = sqrt(dx*sum((U - sol(x')).^2));
p2q = log2(err2_dxq/err2_dxhalfq);