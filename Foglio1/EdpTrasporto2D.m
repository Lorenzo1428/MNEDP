classdef EdpTrasporto2D

    methods
        
        function [U,dt,X,Y] = transp2D(obj,u,X,Y,x,dx,a1,a2,T,l,idProb,idMethod)
            u = obj.bound(u,x,idProb);
            dt = obj.cfl(dx,a1,a2,l,idMethod);
            Nt = floor(T/dt) + 1;
            dt = T/Nt;
            switch idMethod
                case 1
                    U = obj.Upwind2(u,x,dx,dt,Nt,a1,a2,idProb);
                case 2
                    U = obj.LaxF2(u,x,dx,dt,Nt,a1,a2,idProb);
                case 3
                    U = obj.LaxW2(u,x,dx,dt,Nt,a1,a2,idProb);
            end
        end

        function U = Upwind2(obj,u,x,dx,dt,Nt,a1,a2,idProb)
            lambda = dt/dx;
            U =  u;
            for n = 1:Nt
                for i = 2:size(u,1)-1
                    for j = 2:size(u,2)-1
                        U(i,j) = u(i,j) - 0.5*lambda*a1(i,j)*(u(i+1,j) - u(i-1,j)) - 0.5*lambda*a2(i,j)*(u(i,j+1) - u(i,j-1)) + ...
                            0.5*abs(a1(i,j))*lambda*( u(i+1,j) - 2*u(i,j) + u(i-1,j)) + 0.5*abs(a2(i,j))*lambda*( u(i,j+1) - 2*u(i,j) + u(i,j-1) );
                    end
                end
                U = obj.bound(U,x,idProb);
                u=U;
            end
            U = u;
        end

        function U = LaxF2(obj,u,x,dx,dt,Nt,a1,a2,idProb)
            lambda = dt/dx;
            U = u;
            for n = 1:Nt
                for i = 2:size(u,1)-1
                    for j = 2:size(u,2)-1
                        U(i,j) = 0.25*(u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1)) - 0.5*lambda*a1(i,j)*(u(i+1,j) - u(i-1,j)) - 0.5*lambda*a2(i,j)*(u(i,j+1) - u(i,j-1));
                    end
                end
                U = obj.bound(U,x,idProb);
                u = U;
            end   
            U = u;
        end
                
        function U = LaxW2(obj,u,x,dx,dt,Nt,a1,a2,idProb)
            lambda = dt/dx;
            U = u;
            for n = 1:Nt
                for i = 2:size(u,1)-1
                    for j = 2:size(u,2)-1
                        U(i,j) = u(i,j) - 0.5*lambda*a1(i,j)*(u(i+1,j) - u(i-1,j)) - 0.5*lambda*a2(i,j)*(u(i,j+1) - u(i,j-1)) + ...
                            0.5*(a1(i,j)*lambda)^2*( u(i+1,j) - 2*u(i,j) + u(i-1,j)) + 0.5*((a2(i,j))*lambda)^2*( u(i,j+1) - 2*u(i,j) + u(i,j-1) ) + ...
                            0.5*a1(i,j)*a2(i,j)*lambda^2*( u(i+1,j+1) - u(i-1,j+1) - u(i+1,j-1) + u(i-1,j-1));
                    end
                end
                U = obj.bound(U,x,idProb);
                u = U;
            end
            U = u;
        end

    end

    methods(Static)
        function u = bound(u,x,idProb)
            switch idProb 
                case 1
                    u(1,1:end-1) = u(end-1,1:end-1);
                    u(2:end-1,1) = u(2:end-1,end-1);
                    u(end,1:end-1) = u(2,1:end-1);
                    u(1:end,end) = u(1:end,2);
                case 2
                    u(2:end-1,1) = max(1 - 5*x'.^2,zeros(length(x),1));
                    u(2:end-1,end) = 2*u(2:end-1,end-1) - u(2:end-1,end-2);
                    u(1,1:end) = 2*u(2,1:end) - u(3,1:end);
                    u(end,1:end) = 2*u(end-1,1:end) - u(end-2,1:end);
                case 3
                    u(end,2:end-1) = 0;
                    u(2:end,1) = 2*u(2:end,2) - u(2:end,3);
                    u(2:end,end) = 2*u(2:end,end-1) - u(2:end,end-2);
                    u(1,1:end) = 2*u(2,1:end) - u(3,1:end);
            end
        end

        function [u,dx,X,Y,x,a1,a2] = InitCond(dx,idProb)
            
            switch idProb
                case 1
                    N = floor(2*pi/dx) + 1;
                    dx = 2*pi/N;
                    x = -pi:dx:pi;
                    [X,Y] = ndgrid(x,x);
                    u = sin(X).*cos(Y);
                    u = [u,zeros(size(u,1),1)];
                    u = [u; zeros(1,size(u,2))];

                    a1 = [ones(length(x)),zeros(size(X,1),1)];
                    a1 = [a1; zeros(1,size(a1,2))];
                    a2 = [ones(length(x)),zeros(size(X,1),1)];
                    a2 = [a2; zeros(1,size(a1,2))];
                case 2
                    N = floor(4/dx) + 1;
                    dx = 4/N;
                    x = -2:dx:2;
                    [X,Y] = ndgrid(x,x);
                    u = zeros(length(x));
                    u = [u,zeros(size(u,1),1)];
                    u = [u; zeros(1,size(u,2))];
                    u = [zeros(1,size(u,2));u];

                    a1 = [X,zeros(size(X,1),1)];
                    a1 = [a1;zeros(1,size(a1,2))];
                    a1 = [zeros(1,size(a1,2));a1];

                    a2 = [ones(length(x)),zeros(size(X,1),1)];
                    a2 = [a2; zeros(1,size(a2,2))];
                    a2 = [zeros(1,size(a2,2));a2];

                case 3
                    N = floor(4/dx) + 1;
                    dx = 4/N;
                    x = -2:dx:2;
                    [X,Y] = ndgrid(x,x);
                    u = max(0.5 - X.^2 - Y.^2,0);
                    u = [u,zeros(size(u,1),1)];
                    u = [zeros(1,size(u,2));u];
                    u = [zeros(size(u,1),1), u];

                    a1 = [-ones(length(x)),zeros(size(X,1),1)];
                    a1 = [zeros(1,size(a1,2));a1];
                    a1 = [zeros(size(a1,1),1),a1];

                    a2 = [Y,zeros(size(Y,1),1)];
                    a2 = [zeros(1,size(a2,2));a2];
                    a2 = [zeros(size(a2,1),1),a2];
                case 4
                    N = floor(1/dx);
                    dx = 4/N;
                    x = -1:dx:1;
                    y = 0:dx:2;
                    [X,Y] = ndgrid(x,y);
                case 5
            end
        end

        function dt = cfl(dx,a1,a2,l,idMethod)
            switch idMethod
                case 1
                    dt = l*dx/max(max(abs(a1) + abs(a2)));
                case 2
                    dt = l*0.5*dx/max(max(max(abs(a1)),max(max(abs(a2)))));
                case 3
                    dt = l*dx/max(max(abs(a1) + abs(a2)));
            end
        end
    end

end

