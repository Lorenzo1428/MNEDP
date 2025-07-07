clc 
clear
close all

E = EdpTrasporto1D;
T = 2;
a = -1;
vel = @(x) (-0.7*ones(length(x),1));
l = 1;
idMethod = 1;
idInit = 5;

dx = 0.01;
[u0,x,dx,isPer] = E.InitCond(dx,idInit);
a = vel(x);
[U,dt] = E.RigidTrans_d1(dx,T,l,a,u0,idMethod,isPer,idInit);   
Us = sol(x,T,a,idInit);
%err = abs(Us' - U(end,2:end-1));
%err1 = max(err);
f1 = figure;
plot(x,U(end,1:end-1))
ylim([0 1]);
% subplot(1,2,1)
% plot(x,U(end,2:end-1),x,Us);
% legend("aprrox","sol");
% subplot(1,2,2)
% plot(x,err);

% M = 7;
% dx = 0.5;
% for i = 1:M
%     [u0,x,dx,isPer] = E.InitCond(dx,idInit);
%     a = vel(x);
%     [U,dt] = E.RigidTrans_d1(dx,T,l,a,u0,idMethod,isPer,idInit);   
%     Us = sol(x,T,a,idInit);
%     err = abs(Us' - U(end,1:end));
%     errinf(i) = max(err);
%     err1(i) = dx*sum(err);
%     if i ~= 1
%         pinf(i) = log2(errinf(i-1)/errinf(i));
%         p1(i) = log2(err1(i-1)/err1(i));
%     end
%     dx = 0.5*dx;
% end
% f2 = figure;
% semilogy(1:M,errinf,1:M,err1)
% legend("Linf","L1");
% pinf
% p1


% dx1 = 0.5*dx;
% [u0,x,dx1,isPer] = E.InitCond(dx1,idInit);
% [U1,dt1] = E.RigidTrans_d1(dx1,T,l,a,u0,idMethod,isPer,idInit);
% 
% Us = sol(x,T,a,idInit);
% err2 = max(abs(Us' - U1(end,1:end-1)));
% 
% p = log2(err/err2);
% 
% subplot(2,1,1)
% plot(x,U1(end,1:end-1),x,Us)
% legend("aprrox","sol");
% subplot(2,1,2)
% plot(x,abs(U1(end,1:end-1) - Us'))
% 
% t = 0;
% for n = 1:size(U1,1)
%     Us = sol(x,t,a,idInit);
%     plot(x,U1(n,1:end-1),x,Us)
%     legend("aprrox","sol");
%     ylim([0 1])
%     xlim([-2 2])
%     title("t= " + t);
%     t = t+dt1;
%     drawnow
%     pause(0)
% end



function U = sol(x,t,a,idInit)
    switch idInit
        case 1 % I = [-1 3]
            U = (sin(2*pi*(x-a.*t)).*(x - a*t >= -1 & x - a*t <= 1));
        case 2
            U = (sin(4*pi*(x-a.*t)).*(x - a*t >= -1 & x - a*t <= 0));
        case 3 % I = [0 7]
            U = sin(t - (x - 7)./a).*(x - a*t >= 7);
        case 4 % I = [-2 2]
            U = (1 - (x - a*t).^2).*(x - a*t >= -1 & x - a*t <= 1);
        case 5
            U = (1 - abs((x - a.*t))).*(x - a.*t >= -1 & x - a.*t <= 1);
    end
end