function [u,it] = fivePoints(dx,X,Y,N,f,phi)    
    F = f(X,Y);
    u = phi(X,Y);
    it_max = 1e6;
    tol = 1e-12;
    it = 0;
    eps = 1;
    while eps > tol && it < it_max
        prev = u;
        for i = 2:N-1
            for j = 2:N-1
                u(i,j) = 0.25*(dx*dx*F(i,j) + u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1));
            end
        end
        it = it + 1;
        eps = max(max(abs(u - prev)));
    end
end

