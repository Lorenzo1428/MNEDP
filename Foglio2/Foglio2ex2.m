clc
clear
close all

prompt_epsi = "Scegliere eps: " ;
epsi = input(prompt_epsi);

vel1 = @(x,y) 2*y.*(1 - x.^2);
vel2 = @(x,y) -2*x.*(1 - y.^2);

prompt_prob = "1. Condizioni di Dir 2. Condizioni miste 3. Condizioni miste 2 " + newline + "Scegli il problema: ";
id = input(prompt_prob);

switch id
    case 1
        phi = @(x,y) 1.*(y == -1);
        csi = @(x,y) zeros(length(x),length(y));
    case 2 
        phi = @(x,y) 1.*(x == -1);
        csi = @(x,y) zeros(length(x),length(y));
    case 3
        phi = @(x,y) 1.*(x == -1);
        csi = @(x,y) sin(pi*y).*(x == 1);
end

dx = 0.05;
x = -1:dx:1;
N = length(x);
[X,Y] = ndgrid(x);
a1 = vel1(X,Y);
a2 = vel2(X,Y);
neu = csi(X,Y);
U = phi(X,Y);
if id == 1
    U1 = DiffTransDir(dx,N,epsi,a1,a2);
    U1 = reshape(U1,N-2,N-2);
    U(2:end-1,2:end-1) = U1;
else
    U1 = DiffTransMix(dx,N,neu,epsi,a1,a2);
    U1 = reshape(U1,N-2,N-1);
    U(2:end,2:end-1) = U1';
end
surf(X,Y,U)
title("Soluzione approssimata");



function U = DiffTransDir(dx,N,epsi,a1,a2) 
    A = zeros((N-2)*(N-2));
    F = zeros((N-2)*(N-2),1);
    for i = 2:N-1
        D = diag(4*epsi + dx*abs(a1(2:end-1,i)) + dx*abs(a2(2:end-1,i))); % diagonal block
        E1 = diag(-epsi + 0.5*dx*a1(3:end-1,i) - 0.5*dx*abs(a1(3:end-1,i)) ,1);  %sup diag diagonal block
        F1 = diag(-epsi - 0.5*dx*a1(2:end-2,i) - 0.5*dx*abs(a1(2:end-2,i)),-1); %inf diag diagonal block
        j = i-1;
        A((j-1)*(N-2) +1:j*(N-2),(j-1)*(N-2) +1:j*(N-2)) = D + E1 + F1;
        if j < N-2
            E2 = diag(-epsi + 0.5*dx*a2(2:end-1,i+1) - 0.5*dx*abs(a2(2:end-1,i+1)) );
            F2 = diag(-epsi - 0.5*dx*(a2(2:end-1,i)) - 0.5*dx*abs(a2(2:end-1,i)) );
            A((j-1)*(N-2) +1:j*(N-2),j*(N-2) +1:(j+1)*(N-2)) = E2;
            A(j*(N-2) +1:(j+1)*(N-2),(j-1)*(N-2) +1:j*(N-2)) = F2;
        end
    end
    F(1:N-2) = epsi + 0.5*dx*a2(2:N-1,1) + 0.5*dx*abs(a2(2:N-1,1));
    U = A\F;
end

function U = DiffTransMix(dx,N,csi,epsi,a1,a2) 
    A = zeros((N-1)*(N-2));
    F = zeros((N-1)*(N-2),1);
    for i = 2:N-1
        D = diag(4*epsi + dx*abs(a1(i,2:end-1)) + dx*abs(a2(i,2:end-1))); % diagonal block
        E1 = diag(-epsi + 0.5*dx*a2(i,3:end-1) - 0.5*dx*abs(a2(i,3:end-1)) ,1);  %sup diag diagonal block
        F1 = diag(-epsi - 0.5*dx*a2(i,2:end-2) - 0.5*dx*abs(a2(i,2:end-2)),-1); %inf diag diagonal block
        j = i-1;
        A((j-1)*(N-2) +1:j*(N-2),(j-1)*(N-2) +1:j*(N-2)) = D + E1 + F1;
        if j < N-2
            E2 = diag(-epsi + 0.5*dx*a1(i+1,2:end-1) - 0.5*dx*abs(a1(i+1,2:end-1)) );
            F2 = diag(-epsi - 0.5*dx*(a1(i,2:end-1)) - 0.5*dx*abs(a1(i,2:end-1)) );
            A((j-1)*(N-2) +1:j*(N-2),j*(N-2) +1:(j+1)*(N-2)) = E2;
            A(j*(N-2) +1:(j+1)*(N-2),(j-1)*(N-2) +1:j*(N-2)) = F2;
        end
    end
    D = diag(ones(N-2,1)); % neumann block
    j = N-1;
    A((j-1)*(N-2) +1:j*(N-2),(j-1)*(N-2) +1:j*(N-2)) = D;
    j = N-2;
    A(j*(N-2) +1:(j+1)*(N-2),(j-1)*(N-2) +1:j*(N-2)) = -D;
    F(1:N-2) = epsi + 0.5*dx*a1(1,2:N-1) + 0.5*dx*abs(a1(1,2:N-1));
    F(end-N+3:end) = dx*csi(end,2:end-1);
    U = A\F;
end
