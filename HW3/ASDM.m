function [x, x0] = ASDM(f, df, x_initial, x0_intial, a, b,BETA, MAX_ITER, TOL, debug)



    x = x_initial;
    x0 = x0_intial;
    
    x_cat_prev = [x; x0];
    x_hat_prev = x_cat_prev;
    
    lambda_prev = 0;
    
    
    iter = 1
    fvals = [];
    gvals = [];
    fvals(iter) = f(x,x0, a,b);
    gvals(iter) = norm(df(x,x0,a,b));
    
    progress = @(iter,x,x0, fvals, delta_f,delta_x, lambda, alpha) fprintf('\n------------------ASDM iter = %d:------------------\n (x = %s, x_0=%f, F(x)=%f delta_F = %f,delta_x = %f, lambda= %f, alpha = %f) \n -------------------------------------------- \n ' , ...
    iter, mat2str(x,6),num2str(x0), fvals, delta_f, delta_x, lambda, alpha);

       
    if debug
            disp(sprintf('-----------------------Iteration: %d--------------------------------', iter));
            disp('     x_1        x_2        x_0     f(x)     delta_F   delta_x   lambda    alpha')
            disp([  x(1), x(2), x0, fvals(iter), NaN,  NaN, 0, NaN])
            
    end
    
    norm_grad = 10000
    delta_f = 10000;
    delta_x = 10000;
    alpha = 10000;
    
    
    lambda_prev = (1 + (1 + 4*lambda_prev^2)^(.5))/2;
    while iter < MAX_ITER 
        if abs(norm_grad) < TOL
            disp(sprintf('----GRADIENT NORM IS BELOW TOLERANCE CONVERGENCE OF FUNCTION AFTER %d ITERATIONS-----', iter))
            disp(sprintf('---------------------------FINAL Iteration: %d----------------------------------', iter));
            disp('     x_1        x_2        x_0     f(x)     delta_F   delta_x   lambda    alpha   NORMGRAD')
            disp([  x(1), x(2), x0, fvals(iter),  delta_f,  delta_x, lambda_new, alpha, norm_grad])
            break;
        end
        
        if delta_x < 1e-8 
            iter = iter + 1;
            disp(sprintf('CHANGE IN X IS TINY, CONVERGENCE OF FUNCTION AFTER %d ITERATIONS', iter))
            disp(sprintf('---------------------------FINAL Iteration: %d----------------------------------', iter));
            disp('     x_1        x_2        x_0     f(x)     delta_F   delta_x   lambda    alpha   NORMGRAD')
            disp([  x(1), x(2), x0, fvals(iter),  delta_f,  delta_x, lambda_new, alpha, norm_grad])
            break;
        end
        
        iter = iter + 1;
        lambda_new = (1 + (1 + 4*lambda_prev^2)^(.5))/2;
        alpha = (1 - lambda_prev)/(lambda_new);
        disp(alpha)
        
        x_hat_new = x_cat_prev - (1/BETA)*df(x,x0, a, b);
        x_cat_new = (1 - alpha)*x_hat_new + alpha*x_hat_prev;
             
        x = x_cat_new(1:end-1);
        x0 = x_cat_new(end);
        
        fvals(iter) = f(x,x0, a, b); 
        gvals(iter) = norm(df(x,x0,a,b));
        
        delta_f = fvals(iter) - fvals(iter - 1); 
        delta_x = norm(x_cat_new - x_cat_prev, 2);
        norm_grad = norm(df(x,x0,a,b),2); 
        
        

        
        if debug
            disp(sprintf('-----------------------Iteration: %d--------------------------------', iter));
            disp('     x_1        x_2        x_0     f(x)     delta_F   delta_x   lambda    alpha   NORMGRAD')
            disp([  x(1), x(2), x0, fvals(iter),  delta_f,  delta_x, lambda_new, alpha, norm_grad])
        end
    
        x_cat_prev = x_cat_new ;
        x_hat_prev = x_hat_new;
        lambda_prev = lambda_new; 
    end
    
    figure(); 
    plot(1:iter, gvals, 'LineWidth',2); grid on;
    title('Norm Gradient of Objective Function ASDM'); xlabel('Iteration'); ylabel('g(x)');
    
    




end
