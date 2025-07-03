classdef EdpTrasporto2D

    methods
        
        function [U,dt,X,Y] = laxF2(obj,dx,a,T)
            N = floor(2*pi/dx) + 1;
            dx = 2*pi/N;
            x = -pi:dx:pi;
            [u,X,Y] = obj.initCond(x);
            u = obj.bound(u);
            U = u;
            dt = obj.cfl(dx,a);
            Nt = floor(T/dt) + 1;
            dt = T/Nt;
            lambda = dt/dx;
            for n = 1:Nt
                %u(2:end-1,2:end-1) = 0.25*(u(3:end,2:end-1) + u(1:end-2,2:end-1) + u(2:end-1,3:end) + u(2:end-1,1:end-2)) - 0.5*lambda*a(1)*(u(3:end,2:end-1) - u(1:end-2,2:end-1)) - 0.5*lambda*a(2)*(u(2:end-1,3:end) - u(2:end-1,1:end-2));
                for i = 2:size(u,1)-1
                    for j = 2:size(u,2)-1
                        U(i,j) = 0.25*(u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1)) - 0.5*lambda*a(1)*(u(i+1,j) - u(i-1,j)) - 0.5*lambda*a(2)*(u(i,j+1) - u(i,j-1));
                        %U(i,j) = u(i,j) - 0.5*a(1)*lambda*(u(i+1,j) - u(i-1,j)) - 0.5*a(2)*lambda*(u(i,j+1) - u(i,j-1)) + 0.5*(lambda*a(1))^2*(u(i+1,j) - 2*u(i,j) + u(i-1,j)) + 0.5*(lambda*a(2))^2*(u(i,j+1) - 2*u(i,j) + u(i,j-1)) + 0.5*a(1)*a(2)*(lambda)*lambda*(u(i+1,j+1) - u(i-1,j+1) - u(i+1,j-1) + u(i-1,j-1));
                    end
                end
                U = obj.bound(U);
                u = U;
                % surf(X,Y,U(1:end-1,1:end-1));
                % zlim([-1 1])
                % drawnow;
                % pause(0.1)
            end
        end

    end

    methods(Static)
        function u = bound(u)
            u(1,1:end-1) = u(end-1,1:end-1);
            u(2:end-1,1) = u(2:end-1,end-1);
            u(end,2:end-1) = u(2,2:end-1);
            u(2:end-1,end) = u(2:end-1,2);
        end

        function [u,X,Y] = initCond(x)
            [X,Y] = ndgrid(x,x);
            u = sin(X).*cos(Y);
            u = [u,zeros(size(u,1),1)];
            u = [u; zeros(1,size(u,2))];
        end

        function dt = cfl(dx,a)
            dt = 0.8*dx/max(a(1),a(2));
            %dt = 0.1*dx;
        end
    end

end

