function [x_k, iter, resid] = mysd(matA, b, x_0, st_limit, iter_limit)

iter = 0;
resid = [];
x_k = x_0;
normb = norm(b); 
residual = norm(b - matA * x_0);
resid(1) = residual/normb; 
rk =  b - matA * x_k;

  while(residual/normb > st_limit && iter < iter_limit) 
    alpha_k = (rk' * rk)/(rk' * matA * rk);
    x_k = x_k + alpha_k * rk; 
    iter = iter + 1;
    rk =  b - matA * x_k;
    residual = norm(rk);
    resid = [resid residual/normb];
  end

end
    
    
    