clc
clear
close all

dx = 0.1;

a = -pi;

N = floor(2*pi/dx); 
x = a + dx*(0:N)';

[X, Y] = meshgrid(x,x);

z = sin(X).*cos(Y);
surf(X,Y,z)
