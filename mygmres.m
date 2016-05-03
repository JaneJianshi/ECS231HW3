function [x_0, iter, resid] = mygmres(matA, b, x_0, st_limit, iter_limit, restart)

iter = 0;
resid = [];

normb = norm(b); 
n = size(matA, 2);

r0 =  b - matA * x_0;
beta = norm(r0);
resid(1) = beta/normb;
%v_1 = r0 / beta; 


  while(beta/normb > st_limit && iter < iter_limit) 
    m = restart;
    v_1 = r0 / beta;
    V(1:n,1:(m+1)) = zeros(n,m+1);
    V(:, 1) = v_1;
    H(1:(m+1),1:m) = zeros(m+1,m);
    for j = 1: m
        w = matA * V(:, j);
        for i = 1:j
            H(i, j) = V(:, i)' * w; 
            w = w - H(i,j) * V(:, i); 
        end
        
        H(j+1, j) = norm(w); 
        if (H(j+1, j) == 0)
            break; 
        end
        
        V(:, j+1) = w/H(j+1, j);
    end
    
    e1 = zeros(m+1,1);
    e1(1) = 1.0;
    y = H(:, 1:m)\(beta * e1); 
    
    x_0 = x_0 + V(: , 1:m) * y; 
    r0 =  b - matA * x_0;
    
    
    iter = iter + 1;
    
    beta = norm(r0);
    resid = [resid beta/normb];
  end

end