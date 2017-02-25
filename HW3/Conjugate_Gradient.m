function [x, x0] = ASDM(f, df, x_initial, x0_intial, a, b,BETA, MAX_ITER, TOL)



    x = x_initial;
    x0 = x0_intial;
    
    x_cat_prev = [x; x0];
    x_hat_prev = x_cat_prev;
    
    lambda_prev = 0;
    
    
    iter = 1
    fvals = [];
    fvals(iter) = f(x,x0, a,b);
    
    progress = @(iter,x,x0, delta_f,delta_x, lambda, alpha) ... 
    fprintf('------------------ASDM iter = %3d:------------------\n (x = %s, x_0=%f, F(x)=%f delta_F = %f,delta_x = %f, lambda= %f, alpha = %f) \n -------------------------------------------- \n ' , ...
    iter, mat2str(x,6),num2str(x0),  f(x,x0, a, b), delta_f, delta_x, lambda, alpha);

    progress(iter,x,x0, 1000000,100000, lambda_prev, 0);
    delta_f = 1000000000;
    
    lambda_prev = (1 + (1 + 4*lambda_prev^2)^(.5))/2;
    while iter < MAX_ITER && abs(delta_f) > TOL
        iter = iter + 1;
        lambda_new = (1 + (1 + 4*lambda_prev^2)^(.5))/2;
        alpha = (1 - lambda_prev)/(lambda_new);
        disp(alpha)
        
        x_hat_new = x_cat_prev - (1/BETA)*df(x,x0, a, b);
        x_cat_new = (1 - alpha)*x_hat_new + alpha*x_hat_prev;
             
        x = x_cat_new(1:end-1);
        x0 = x_cat_new(end);
        
        fvals(iter) = f(x,x0, a, b); 
        
        delta_f = fvals(iter) - fvals(iter - 1); 
        delta_x = norm(x_cat_new - x_cat_prev, 2);
  
        progress(iter, x,x0, delta_f, delta_x, lambda_new, alpha); 
        x_cat_prev = x_cat_new ;
        x_hat_prev = x_hat_new;
        lambda_prev = lambda_new; 
    end
    
    figure(); 
    plot(1:iter, fvals, 'LineWidth',2); grid on;
    title('Objective Function ASDM'); xlabel('Iteration'); ylabel('F(x)');
    
    




end
