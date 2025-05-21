clc 
clear
close all

T = EdpTrasportoD1;
T.time = 0.6;
T.velocity = 3;
%T.velocity = @(x) (-x);

a = T.velocity;
T.amplitude = 0.5;

dx = 0.01;
[U,u0,x,dt,tf] = T.RigidTrans_d1(dx,1,4);

Us = sol(x - a*(tf));
err = norm(Us' - U(end,1:end-1),inf);


dx1 = 0.5*dx;
[U1,u0,x,dt1,tf1] = T.RigidTrans_d1(dx1,1,4);

Us = sol(x - a*(tf1));
err2 = norm(Us' - U1(end,1:end-1),inf);

log(err/err2)

plot(x,U1(end,1:end-1),x,Us)

% t = 0;
% for n = 1:size(U,1)
%     Us = sol(x - a*t);
%     plot(x,U(n,1:end-1),x,Us)
%     ylim([-2 2])
%     xlim([-1 3])
%     t = t+dt;
%     title("t= " + t)
%     drawnow
%     pause(0)
% end



function U = sol(x)
    U = (1 - x.^2).*(x >= -1 & x <= 1);
end