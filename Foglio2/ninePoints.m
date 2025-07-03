function [u,it] = ninePoints(dx,X,Y,N,f,phi)    
    F = f(X,Y);
    u = phi(X,Y);
    it_max = 1e6;
    tol = 1e-20;
    it = 0;
    eps = 1;
    while eps > tol && it < it_max
        prev = u;
        for i = 2:N-1
            for j = 2:N-1
                u(i,j) = dx*dx*(1/40)*( F(i+1,j) + F(i-1,j) + F(i,j+1) + F(i,j-1) + 8*F(i,j) ) + ... 
                    (1/20)*( u(i+1,j+1) + u(i-1,j+1) + u(i+1,j-1) + u(i-1,j-1) ) +...
                    (1/5)*( u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1));
            end
        end
        it = it + 1;
        eps = max(max(abs(u - prev)));
    end
end
