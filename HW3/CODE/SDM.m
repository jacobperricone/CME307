%%
% Author: Jacob Perricone
% Description: Implements the Steepest Descent Method on an initial point x0_initial, x_initial
% for classifying two sets of points a and b
% f: function to optimize over
% df: gradient function
%%


function [x, x0, iter,fvals, gvals] = SDM(f, df, x_initial, x0_intial, a, b,ALPHA, MAX_ITER, TOL, debug,mu)



x = x_initial;
x0 = x0_intial;

x_cat_prev = [x; x0];

iter = 1;
fvals = [];
gvals = [];
gvals(iter) = norm(df(x,x0,a,b,mu));
fvals(iter) = f(x,x0, a,b,mu);


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
    x_cat_new = x_cat_prev - ALPHA*df(x,x0, a, b,mu);
    
    x = x_cat_new(1:end-1);
    x0 = x_cat_new(end);
    
    fvals(iter) = f(x,x0, a, b,mu);
    gvals(iter) = norm(df(x,x0,a,b,mu));
    
    delta_f = fvals(iter) - fvals(iter - 1);
    delta_x = norm(x_cat_new - x_cat_prev, 2);
    norm_grad = norm(df(x,x0,a,b,mu),2);
    
    
    if debug
        disp(sprintf('-----------------------Iteration: %d--------------------------------', iter));
        disp('     x_1        x_2        x_0     f(x)     delta_F   delta_x   NORMGRAD')
        disp([ x(1), x(2), x0, fvals(iter), delta_f,  delta_x, norm_grad])
        
    end
    x_cat_prev = x_cat_new ;
end

if mu
    s = 'With_Regularization';
else
    s = 'No_Regularization';
end


fig = figure();
subplot(2,1,1)
plot(1:iter, gvals, 'LineWidth',2); grid on;
title(strcat({' Norm Gradient Function SDM Method: '}, s));
xlabel('Iteration'); ylabel('G(x)');

subplot(2, 1, 2)
plot(1:iter, fvals, 'LineWidth',2); grid on;
title(strcat({'Objective Function SDM Method: '}, s));

xlabel('Iteration'); ylabel('F(x)');

save_string = [pwd strcat('/FIGURES/Problem5/SDM', '_', s,'.png')]
saveas(fig, save_string)



end
