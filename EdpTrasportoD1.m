classdef EdpTrasportoD1

    properties
        velocity = 1;
        time = 10;
        amplitude = 0.1;
        timeStep = 0.1;
    end

    methods

        function [u0,xn,isPer] = InitCond(obj,dx,idInit)
            
            switch idInit
                case 1 % I = [-1 3]
                    N = floor(4/dx);
                    xn = -1 + (0:N)'.*dx;
                    u0 = (sin(2*pi*xn).*(xn >= -1 & xn <= 1));
                    isPer = true;
                case 2
                    N = floor(4/dx);
                    xn = -1 + (0:N)'.*dx;
                    u0 = (sin(4*pi*xn).*(xn >= -1 & xn <= 0));
                    isPer = true;
                case 3 % I = [0 7]
                    N = floor(7/dx);
                    xn = (0:N)'.*dx;
                    u0 = zeros(N+1,1);
                    isPer = false;
                case 4 % I = [-2 2]
                    N = floor(4/dx);
                    xn = -2 + (0:N)'.*dx;
                    u0 = (1 - xn.^2).*(xn >= -1 & xn <= 1);
                    isPer = false;
                case 5
                    N = floor(4/dx);
                    xn = -2 + (0:N)'.*dx;
                    u0 = (1 - abs(xn)).*(xn >= -1 & xn <= 1);
                    isPer = false;
            end
        end
        
        function [U,u0,nodes,dt,t] = RigidTrans_d1(obj,dx,idMethod,idInit)
            T = obj.time;
            a = obj.velocity;
            %dx = obj.amplitude;
            [u0,nodes,isPer] = obj.InitCond(dx,idInit);

            u = obj.SetBounds(u0,a);
            u = obj.BoundCond(u,0,a,idInit,isPer);
            dt = obj.Cfl(dx,a); 
            
            U(1,:) = u;
            
            switch idMethod
                case 1
                    [U,t] = obj.Upwind1(u,U,T,a,a,dt,dx,idInit,isPer);
                case 2
                    disp(2)
                    [U,t] = obj.LaxF1(u,U,T,a,a,dt,dx,idInit,isPer);
                case 3
                    [U,t] = obj.LaxW1(u,U,T,a,a,dt,dx,idInit,isPer);
            end
                
        end

        function [U,u0,nodes,dt] = NonRigidTrans_d1(obj,idMethod,idInit)
            T = obj.time;
            a = obj.velocity;
            dx = obj.amplitude;
            [u0,nodes,isPer] = obj.InitCond(idInit);
            
            a = a(nodes);
            [u,a1] = obj.SetBounds(u0,a);
            u = obj.BoundCond(u,0,a,idInit,isPer);
            dt = obj.Cfl(dx,a); 
            
            U(1,:) = u';
            
            switch idMethod
                case 1
                    U = obj.Upwind1(u,U,T,a,a1,dt,dx,idInit,isPer);
                case 2
                    U = obj.LaxF1(u,U,T,a,a1,dt,dx,idInit,isPer);
                case 3
                    U = obj.LaxW1(u,U,T,a,a1,dt,dx,idInit,isPer);
            end
                
        end

        function [U,t] = Upwind1(obj,u,U,T,a,a1,dt,dx,idInit,isPer)
            lambda = dt/dx;
            t = 0;
            n = 2;
            while t <= T    
                u(2:end-1) = u(2:end-1) - 0.5*lambda*a1.*(u(3:end) - u(1:end-2)) + 0.5*lambda*abs(a1).*(u(3:end) -2*u(2:end-1) + u(1:end-2));
                u = obj.BoundCond(u,t,a,idInit,isPer);
                U(n,:) = u;
                n = n+1;
                t = t + dt;
            end
        end

        function [U,t] = LaxF1(obj,u,U,T,a,a1,dt,dx,idInit,isPer)
            t = 0;
            n = 2;
            lambda = (0.5*a1.*dt)/dx;

            while t <= T    
                u(2:end-1) = 0.5*(u(3:end) + u(1:end-2)) - lambda*(u(3:end) - u(1:end-2));
                u = obj.BoundCond(u,t,a,idInit,isPer);
                U(n,:) = u;
                n = n+1;
                t = t + dt;
            end
        end

        function [U,t] = LaxW1(obj,u,U,T,a,a1,dt,dx,idInit,isPer)
            lambda = dt/dx;
            t = 0;
            n = 2;
            while t <= T    
                u(2:end-1) = u(2:end-1) - 0.5*a1.*lambda*(u(3:end) - u(1:end-2)) + 0.5*(a1.*lambda)^2*(u(3:end) -2*u(2:end-1) + u(1:end-2));
                u = obj.BoundCond(u,t,a,idInit,isPer);
                U(n,:) = u;
                n = n+1;
                t = t + dt;
            end
        end
    end

    methods(Static)

        function dt = Cfl(dx,a)
            dt = 0.6*dx/norm(a,inf);
        end

        function u = BoundCond(u,t,a,idInit,isPer)
            if isPer
                u(1) = u(end-1);
                u(end) = u(2);
            elseif ~isPer && idInit == 3
                u(end) = sin(t);
                u(1) = 2*u(2) - u(3);
            else 
                if isscalar(a)
                    if a >= 0
                        u(1) = 0;
                        u(end) = 2*u(end-1) - u(end-2);
                    else
                        u(end) = 0;
                        u(1) = 2*u(2) - u(3);
                    end
                else
                    if a(1) >= 0
                        u(1) = 0;
                    else
                        u(1) = 2*u(2) - u(3);
                    end

                    if a(end) >= 0 
                        u(end) = 0;
                    else
                        u(end) = 2*u(end-1) - u(end-2);
                    end
                end
            end
        end

        function [u,a1] = SetBounds(u,a)
            if isscalar(a)
                if a >= 0
                    u = [u;0];
                else
                    u = [0;u];
                end
            else
                if a(1) >= 0  && a(end) < 0
                    a1 = a(2:end-1);
                elseif a(1) >= 0 && a(end) >=0
                    u = [u;0];
                    a1 = a(2:end);
                elseif a(1) < 0 && a(end) >=0
                    a1 = a;
                    u = [0;u;0];
                else
                    u = [0;u];
                    a1 = a(1:end-1);
                end              
            end
        end

    end

end
