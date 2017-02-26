%%
% Author: Jacob Perricone
% Description: Implements the Conjugate Gradient Method on an initial point x0_initial, x_initial
% for classifying two sets of points a and b
% f: function to optimize over
% df: gradient function
%%
function [x, x0] = CGD(f, df, x_initial, x0_intial, a, b, MAX_ITER, TOL, debug,mu)



% save initial coordinates
x = x_initial;
x0 = x0_intial;


x_cat_prev = [x;x0];
g_prev = df(x,x0,a,b, mu);
d_prev = -g_prev;

iter = 1;
fvals = [];
gvals = [];


fvals(iter) = f(x,x0, a,b,mu);
gvals(iter) = norm(df(x,x0,a,b, mu));
norm_grad = 10000

progress = @(iter,x,x0, delta_f,delta_x, i) fprintf('------------------CGD iter = 6%d:------------------\n (x = %s, x_0=%s, F(x)=%f delta_F = %f,delta_x = %f, i= %f) \n -------------------------------------------- \n ' , ...
    iter, mat2str(x,2),num2str(x0),  f(x,x0, a, b), delta_f, delta_x, i);


if debug
    disp(sprintf('-----------------------Iteration: %d--------------------------------', iter));
    disp('     x_1        x_2        x_0     f(x)     delta_F   delta_x,  gradF i ')
    disp([  x(1), x(2), x0, fvals(iter), NaN,  NaN, 1])
    
end
i = 1;


while i < 3
    
    
    alpha = (-g_prev'*d_prev)/(d_prev'*fvals(iter)*d_prev);
    x_cat_new = x_cat_prev + alpha*d_prev;
    x_cat_prev;
    
    x = x_cat_new(1:end-1);
    x0 = x_cat_new(end);
    
    fvals(iter + 1) = f(x,x0,a,b,mu);
    gvals(iter + 1) = norm(df(x,x0,a,b,mu));
    
    delta_f = fvals(iter+1) - fvals(iter);
    norm_grad =  norm(df(x,x0,a,b,mu),2);
    delta_x = norm(x_cat_new - x_cat_prev, 2);
    
    
    
    if debug
        disp(sprintf('-----------------------Iteration: %d--------------------------------', iter));
        disp('     x_1        x_2        x_0     f(x)     delta_F   delta_x   NORMGRAD  i ')
        disp([  x(1), x(2), x0, fvals(iter + 1), delta_f,  delta_x, norm_grad, i])
        
    end
    
    if abs(norm_grad) < TOL
        iter = iter + 1;
        disp(sprintf('----GRADIENT NORM IS BELOW TOLERANCE CONVERGENCE OF FUNCTION AFTER %d ITERATIONS-----', iter))
        disp(sprintf('-----------------------FINAL RESULT ITERATIONA: %d--------------------------------', iter));
        disp('     x_1        x_2        x_0     f(x)     delta_F   delta_x   NORMGRAD  i ')
        disp([  x(1), x(2), x0, fvals(iter), delta_f,  delta_x, norm_grad, i])
        break;
    end
    
%     if delta_x < 1e-8
%         iter = iter + 1;
%         disp(sprintf('CHANGE IN X IS TINY, CONVERGENCE AFTER %d ITERATIONS', iter))
%         disp(sprintf('-----------------------FINAL RESULT ITERATION: %d--------------------------------', iter));
%         disp('     x_1        x_2        x_0     f(x)     delta_F   delta_x   NORMGRAD  i ')
%         disp([  x(1), x(2), x0, fvals(iter), delta_f,  delta_x, norm_grad, i])
%         break;
%     end
    
    if iter+1 >  MAX_ITER
        iter = iter + 1;
        disp('MAXIMIUM ITERATIONS REACHED \n ')
        break;
    end
    
    
    if i ~= 2
        g_prev;
        
        g_new = df(x,x0,a,b,mu);
        beta = (g_new'*fvals(iter)*d_prev)/(d_prev'*fvals(iter)*d_prev);
        beta ;
        d_prev;
        
        d_new = -g_new + beta*d_prev;
        d_prev = d_new;
        g_prev = g_new;
        x_cat_prev = x_cat_new;
        i = i+1;
        
    else
        if debug
            disp('------------------------------RESETTING----------------------------')
        end
        x_cat_prev = x_cat_new;
        g_new = df(x,x0,a,b,mu);
        g_prev = g_new;
        d_prev = -g_new;
        
        i=1;
        
        
    end
    
    iter = iter + 1;
    
end

if mu
    s = 'With Regularization';
else
    s = 'No Regularization';
end


fig = figure();
subplot(2,1,1)
plot(1:iter, gvals, 'LineWidth',2); grid on;
title(strcat({' Norm Gradient Function CGD Method: '}, s));
xlabel('Iteration'); ylabel('G(x)');

subplot(2, 1, 2)
plot(1:iter, fvals, 'LineWidth',2); grid on;
title(strcat({'Objective Function CGD Method: '}, s));

xlabel('Iteration'); ylabel('F(x)');

save_string = [pwd strcat('/FIGURES/Problem5/CGD', '_', s,'.png')]
saveas(fig, save_string)





end
