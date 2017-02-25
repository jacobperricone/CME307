function [x, x0] = SDM(f, df, x_initial, x0_intial, a, b,ALPHA, MAX_ITER, TOL, debug)



x = x_initial;
x0 = x0_intial;

x_cat_prev = [x; x0];

iter = 1;
fvals = [];
fvals(iter) = f(x,x0, a,b);
gvals = [];
gvals(iter) = norm(df(x,x0,a,b));



delta_f = 1000;
delta_x = 1000;
norm_grad = 10000

if debug
    disp(sprintf('-----------------------Iteration: %d--------------------------------', iter));
    disp('     x_1        x_2        x_0     f(x)     delta_F   delta_x')
    disp([  x(1), x(2), x0, fvals(iter), NaN,  NaN])
    
end


while iter < MAX_ITER
    if abs(norm_grad) < TOL
        disp(sprintf('----GRADIENT NORM IS BELOW TOLERANCE CONVERGENCE OF FUNCTION AFTER %d ITERATIONS-----', iter))
        disp(sprintf('-----------------------FINAL Iteration: %d--------------------------------', iter));
        disp('     x_1        x_2        x_0     f(x)     delta_F   delta_x   NORMGRAD')
        disp([ x(1), x(2), x0, fvals(iter), delta_f,  delta_x, norm_grad])
        break;
    end
    
    if delta_x < 1e-8
     
        disp(sprintf('CHANGE IN X IS TINY, CONVERGENCE OF FUNCTION AFTER %d ITERATIONS', iter))
        disp(sprintf('-----------------------FINAL Iteration: %d--------------------------------', iter));
        disp('     x_1        x_2        x_0     f(x)     delta_F   delta_x   NORMGRAD')
        disp([ x(1), x(2), x0, fvals(iter), delta_f,  delta_x, norm_grad])
        break;
    end
    
    iter = iter + 1;
    x_cat_new = x_cat_prev - ALPHA*df(x,x0, a, b);
    
    x = x_cat_new(1:end-1);
    
    x0 = x_cat_new(end);
    
    fvals(iter) = f(x,x0, a, b);
    gvals(iter) = norm(df(x,x0,a,b));
    
    delta_f = fvals(iter) - fvals(iter - 1);
    delta_x = norm(x_cat_new - x_cat_prev, 2);
    norm_grad = norm(df(x,x0,a,b),2);
    
    
    if debug
        disp(sprintf('-----------------------Iteration: %d--------------------------------', iter));
        disp('     x_1        x_2        x_0     f(x)     delta_F   delta_x   NORMGRAD')
        disp([ x(1), x(2), x0, fvals(iter), delta_f,  delta_x, norm_grad])
        
    end
    x_cat_prev = x_cat_new ;
end

figure();
plot(1:iter, fvals, 'LineWidth',2); grid on;
title('Objective Function SDM'); xlabel('Iteration'); ylabel('F(x)');






end
