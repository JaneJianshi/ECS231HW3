function [xk, iter, resid] = mycg(matA, b, x_0, st_limit, iter_limit)

 normb = norm(b);
 xk =  x_0;
 rk = b - matA * xk;
 pk = rk;
 iter = 0;
 resid = [];
 resid(1) = norm(rk)/normb;
 
  while(norm(rk)/normb > st_limit && iter < iter_limit)
   Apk = matA * pk;
   ak = (rk' * rk)/(pk' * Apk);
   xk_new = xk + ak * pk;
   rk_new = rk - ak * Apk;
   
   betak = (rk_new' * rk_new) / (rk' * rk);
   pk_new = rk_new + betak * pk;
   
   xk = xk_new;
   rk = rk_new;
   pk = pk_new;
   iter = iter + 1;
   resid = [resid norm(rk)/normb]; 
   
  end
  
end
  
  