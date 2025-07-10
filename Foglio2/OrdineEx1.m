clc
clear
close all

a = 0;
b = 1;
f = @(x,y) (32*pi*pi*cos(4*pi*x).*sin(4*pi*y));
phi = @(x,y) (sin(4*pi*y));
sol = @(x,y) (cos(4*pi*x).*sin(4*pi*y));

M = 5;
dx = zeros(M,1);
p5 = zeros(M,1);
p9 = zeros(M,1);
err5 = zeros(M,1);
err9 = zeros(M,1);
dx(1) = 0.1;

for i = 1:M
    N = floor((b-a)/dx(i));
    dx1 = (b-a)/N;
    x = a:dx1:b;
    N = length(x);
    [X,Y] = ndgrid(x);
    S = sol(X,Y);
    [u5,~] = fivePoints(dx1,X,Y,N,f,phi);
    [u9,~] = ninePoints(dx1,X,Y,N,f,phi);
    err5(i) = max(max(abs(u5 - S)));
    err9(i) = max(max(abs(u9 - S)));
    if i>1
        p5(i) = log2(err5(i-1)/err5(i));
        p9(i) = log2(err9(i-1)/err9(i));
    end
    if i < M
        dx(i+1) = 0.5*dx(i);
    end
end
table(dx,err5,p5,err9,p9)