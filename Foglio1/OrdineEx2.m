clc
clear
close all

E = EdpTrasporto2D;
sol = @(x,y,a1,a2,T) sin(x-T*a1).*cos(y-T*a2);
M = 5;
idProb = 1;
T = 1; 
l=1;

dx = zeros(M,1);
err = zeros(M,1);
p = zeros(M,1);
dx(1) = 0.1;
for i = 1:3
    idMethod = i;
    if i == 3
        dx(1) = 0.2; 
        M = 4;
        T = 1;
        l = 0.001;
    end
    for j = 1:M
        [u,dx1,X,Y,x,a1,a2] = E.InitCond(dx(j),idProb);      
        [U,dt,X,Y] = E.transp2D(u,X,Y,x,dx1,a1,a2,T,l,idProb,idMethod);
        S = sol(X,Y,a1(1:end-1,1:end-1),a2(1:end-1,1:end-1),T);
        err(j) = max(max(abs(U(1:end-1,1:end-1) - S)));
        if j>1
        p(j) = log2(err(j-1)/err(j));
        end
        if j < M
            dx(j+1) = 0.5*dx(j);
        end
    end    
    switch i
        case 1
            T_UP = table(dx,err,p)
        case 2
            T_LF = table(dx,err,p)
        case 3
            T_LW = table(dx(1:M),err(1:M),p(1:M))
    end
    dx = zeros(M,1);
    dx(1) = 0.1;
end