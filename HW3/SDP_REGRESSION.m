function [X_final, X] = SDP_REGRESSION(f, df, A,X0,mu,x_true,TOL, MAX_ITER, ALPHA,a, debug)


X_prev = X0;
iter = 1;
fvals = [];
fvals(iter) = f(A,X0, mu);
x_init = X_prev(end,1:2);
delta_f = 1000;
delta_x = 1000;
norm_grad = 10000;

if debug
    disp(sprintf('-----------------------Iteration: %d--------------------------------', iter));
    disp('     x_1        x_2        Z     f(x)     delta_F   delta_x   NORMGRAD   x_true_1  x_true_2 ')
    disp([  X_prev(end,1), X_prev(end,2),X_prev(end,end), fvals(iter), NaN,  NaN, NaN, x_true(1), x_true(2)])
    
end


while iter < MAX_ITER
    if abs(norm_grad) < TOL
        disp(sprintf('----GRADIENT NORM IS BELOW TOLERANCE CONVERGENCE OF FUNCTION AFTER %d ITERATIONS-----', iter))
        disp(sprintf('-----------------------FINAL Iteration: %d--------------------------------', iter));
        disp('     x_1        x_2        Z     f(x)     delta_F   delta_x   NORMGRAD   x_true_1  x_true_2 ')
        disp([  X_prev(end,1), X_prev(end,2),X_prev(end,end), fvals(iter), delta_f, delta_x, norm_grad, x_true(1), x_true(2)])
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
    
    
    
    norm_grad = norm(df(A,X_new, mu),2);
    delta_f = fvals(iter) - fvals(iter - 1);
    
    delta_x = norm(X_new - X_prev,2);
    
   
    if debug
        disp(sprintf('-----------------------Iteration: %d--------------------------------', iter));
        disp('     x_1        x_2        Z     f(x)     delta_F   delta_x   NORMGRAD   x_true_1  x_true_2 ')
        disp([  X_new(end,1), X_new(end,2),X_new(end,end), fvals(iter), delta_f, delta_x, norm_grad, x_true(1), x_true(2)])
        
        
        
    end
    
    X_prev = X_new;
end

X_final = X_prev(end,1:2);
X = X_prev



figure();
subplot(2,1,1)
plot(1:iter, fvals, 'LineWidth',2); grid on;
title(strcat({'Objective Function when mu = '}, num2str(mu)) ); xlabel('Iteration'); ylabel('F(x)');
subplot(2,1,2)
tmp = a;
tmp(:,4) = tmp(:,1);
plot(tmp(1,:), tmp(2,:))
hold on
plot(x_true(1), x_true(2),'-*b','MarkerSize',10)
plot(X_final(1),X_final(2),'-or','MarkerSize',10)
plot(x_init(1),x_init(2),'-*g','MarkerSize',10)
title(strcat({'Estimated Location when mu = '}, num2str(mu)))
xlabel('X_1')
ylabel('X_2')
legend('Convex Hull','True Point', 'Estimate', 'Initial Point')
hold off




    
end