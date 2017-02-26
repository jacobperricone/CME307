function [x1_final,x2_final, X] = SDP_REGRESSION_HW2(f, df, A,X0,mu,x_true,TOL, MAX_ITER, ALPHA,a, debug, type)


X_prev = X0;
iter = 1;
fvals = [];
gvals = [];
fvals(iter) = f(A,X0, mu);
gvals(iter) = norm(df(A,X0,mu));
x_init2 = X_prev(end,1:2);
x_init1 = X_prev(end-1,1:2);
delta_f = 1000;
delta_x = 1000;
norm_grad = 10000;
 
if debug
    disp(sprintf('-----------------------Iteration: %d--------------------------------', iter));
    disp('     x1_1        x1_2     x2_1     x2_2      Z     f(x)     delta_F   delta_x   NORMGRAD   x2_true_1  x1_true_2 x2_true_1  x2_true_2')
    disp([  X_prev(end-1,1) X_prev(end-1, 2)  X_prev(end,1), X_prev(end,2), X_prev(end,end), fvals(iter), NaN,  NaN, NaN, x_true(1,1), x_true(2,1)  x_true(1,2) x_true(2,2)])
    
end


while iter < MAX_ITER
    if abs(norm_grad) < TOL
        disp(sprintf('----GRADIENT NORM IS BELOW TOLERANCE CONVERGENCE OF FUNCTION AFTER %d ITERATIONS-----', iter))
        disp(sprintf('-----------------------FINAL Iteration: %d--------------------------------', iter));
        disp('     x1_1        x1_2     x2_1     x2_2      Z     f(x)     delta_F   delta_x   NORMGRAD   x2_true_1  x1_true_2 x2_true_1  x2_true_2')
        disp([  X_prev(end-1,1) X_prev(end-1, 2)  X_prev(end,1), X_prev(end,2),X_prev(end,end), fvals(iter), delta_f, delta_x, norm_grad, x_true(1,1), x_true(2,1)  x_true(1,2) x_true(2,2)])
        break;
    end
    
%     if abs(delta_x) < 1e-8
%          disp(sprintf('CHANGE IN X IS TINY, CONVERGENCE OF FUNCTION AFTER %d ITERATIONS', iter))
%          disp(sprintf('-----------------------FINAL Iteration: %d--------------------------------', iter));
%          disp('     x_1        x_2        Z     f(x)     delta_F   delta_x   NORMGRAD   x_true_1  x_true_2 ')
%          disp([  X_prev(end,1), X_prev(end,2),X_prev(end,end), fvals(iter), delta_f, delta_x, norm_grad, x_true(1), x_true(2)])
%          break
%     end
    
    iter = iter +1;
    alpha = 1/(ALPHA*norm(X_prev, inf)^2);
    
    
    X_new = X_prev - alpha*df(A,X_prev,mu);
    
    fvals(iter) = f(A,X_new, mu);
    gvals(iter) = norm(df(A,X_new,mu));
    
    
    norm_grad = norm(df(A,X_new, mu),2);
    delta_f = fvals(iter) - fvals(iter - 1);
    
    delta_x = norm(X_new - X_prev,2);
    
   
    if debug
        disp(sprintf('-----------------------Iteration: %d--------------------------------', iter));
        disp('     x1_1        x1_2     x2_1     x2_2      Z     f(x)     delta_F   delta_x   NORMGRAD   x2_true_1  x1_true_2 x2_true_1  x2_true_2')
        disp([  X_prev(end-1,1) X_prev(end-1, 2)  X_prev(end,1), X_prev(end,2),X_prev(end,end), fvals(iter), delta_f, delta_x, norm_grad, x_true(1,1), x_true(2,1)  x_true(1,2) x_true(2,2)])
        
        
        
    end
    
    X_prev = X_new;
end

x1_final = X_prev(end-1,1:2);
x2_final = X_prev(end,1:2);





figure();
subplot(3,1,1)
plot(1:iter, fvals, 'LineWidth',2); grid on;
title(strcat({'Objective Function when mu = '}, num2str(mu),{'  (descent = '}, type,')' )) ; 
xlabel('Iteration'); ylabel('F(x)');

subplot(3,1,2)
plot(1:iter, fvals, 'LineWidth',2); grid on;
title(strcat({'Gradient Function when mu = '}, num2str(mu) ,{'  (descent = '}, type,')' )); 
xlabel('Iteration'); ylabel('F(x)');


subplot(3,1,3)
tmp = a;
tmp(:,4) = tmp(:,1);
plot(tmp(1,:), tmp(2,:))
hold on
plot(x_true(1,1), x_true(2,1),'-*r','MarkerSize',10)
plot(x_true(1,2), x_true(2,2),'-*g','MarkerSize',10)
plot(x1_final(1),x1_final(2),'-or','MarkerSize',10)
plot(x2_final(1),x2_final(2),'-og','MarkerSize',10)
plot(x_init1(1),x_init1(2),'-xr','MarkerSize',10)
plot(x_init2(1),x_init2(2),'-xg','MarkerSize',10)
title(strcat({'Estimated Location when mu = '}, num2str(mu)),{' (descent = '}, type,')' )
xlabel('X_1')
ylabel('X_2')
legend('Convex Hull','X_1 True Point', 'X_2 True Point', 'X_1 Estimate', 'X_2 Estimate', 'X_1 Initial Point', 'X_2 Initial Point')
hold off




    
end