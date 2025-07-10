clc
clear
close all
E = EdpTrasporto2D;

idMethod = 3;
idProb = 4;

dx = 0.005;
T = 1; 
l = 0.3;

[u,dx,X,Y,x,a1,a2] = E.InitCond(dx,idProb);
[U,dt,X,Y] = E.transp2D(u,X,Y,x,dx,a1,a2,T,l,idProb,idMethod);
h = surf(X,Y,U(2:end-1,2:end-1));
set(h,'edgecolor','none');
colormap jet
zlim([0,1]);

 
% sol = sin(X-T*a1(1:end-1,1:end-1)).*cos(Y-T*a2(1:end-1,1:end-1));
% f2 = figure;
% surf(X,Y,sol)

% err = abs(U(1:end-1,1:end-1) - sol);
% f3 = figure; 
% mesh(err);

% %errinf = max(err,[],"all");
% errinf = sqrt(dx*dx*(sum(err.^2,"all")));
% 
% dx = 0.5*dx;
% [u,dx,X,Y,x,a1,a2] = E.InitCond(dx,idProb);
% 
% [U,dt,X,Y] = E.transp2D(u,X,Y,x,dx,a1,a2,T,l,idProb,idMethod);
% sol = sin(X-T*a1(1:end-1,1:end-1)).*cos(Y-T*a2(1:end-1,1:end-1));
% err = abs(U(1:end-1,1:end-1) - sol);
% %errinf1 = max(err,[],"all");
% errinf1 = sqrt(dx*dx*(sum(err.^2,"all")));
% 
% p = log2(errinf/errinf1)