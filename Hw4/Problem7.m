clear all
close all

%% Problem 7
% Part (a) and (b)
A = [1, 1, 1; 1, 1, 2; 1, 2, 2];
beta = [.1, 1, 10];
syms x_1 x_2 x_3 fx_1 fx_2 fx_3

for i=1:3
x1(1) = .01;
x2(1) = -.01;
x3(1) = .01;
y(:, 1) = [1;1;1];
fx1(1) = .01;
fx2(1) = -.01;
fx3(1) = .01;
fy(:, 1) = [1;1;1];
    for k=1:99
        x1(k+1) = solve(A(:, 1)'*y(:, k) + beta(i) * A(:, 1)'*(A(:, 1)*x_1 + ...
            A(:, 2)*x2(k) + A(:, 3)*x3(k)) == 0, x_1);

        x2(k+1) = solve(A(:, 2)'*y(:, k) + beta(i) * A(:, 2)'*(A(:, 1)*x1(k+1) + ...
            A(:, 2)*x_2 + A(:, 3)*x3(k)) == 0, x_2);

        x3(k+1) = solve(A(:, 3)'*y(:, k) + beta(i) * A(:, 3)'*(A(:, 1)*x1(k+1) + ...
            A(:, 2)*x2(k+1) + A(:, 3)*x_3) == 0, x_3);

        y(:, k+1)=y(:, k)+beta(i) * (A(:, 1)*x1(k+1) + A(:, 2)*x2(k+1) + A(:, 3)*x3(k+1));

       
        fx1(k+1) = solve(fx_1+A(:, 1)'*fy(:, k) + beta(i) * A(:, 1)'*(A(:, 1)*fx_1 + ...
            A(:, 2)*fx2(k) + A(:, 3)*fx3(k)) == 0, fx_1);

        fx2(k+1) = solve(fx_2+A(:, 2)'*fy(:, k) + beta(i) * A(:, 2)'*(A(:, 1)*fx1(k+1) + ...
            A(:, 2)*fx_2 + A(:, 3)*fx3(k)) == 0, fx_2);

        fx3(k+1) = solve(fx_3+A(:, 3)'*fy(:, k) + beta(i) * A(:, 3)'*(A(:, 1)*fx1(k+1) + ...
            A(:, 2)*fx2(k+1) + A(:, 3)*fx_3) == 0, fx_3);

        fy(:, k+1)=fy(:, k)+beta(i) * (A(:, 1)*fx1(k+1) + A(:, 2)*fx2(k+1) + A(:, 3)*fx3(k+1));
        
        if mod(k,10) == 0
            disp(sprintf('----- beta: %d  ------', beta(i)));
            disp(sprintf('------- Iteration: %d ----------', k));
            disp('      x1        x2        x3')
            disp([  x1(k), x2(k), x3(k)])
            
            disp(sprintf('----- beta: %d  ------', beta(i)));
            disp(sprintf('------- Iteration: %d ----------', k));
            disp('  x1 w/ f    x2 w/ f    x3 w/ f')
            disp([  fx1(k), fx2(k), fx3(k)])
        end

    end
    
    idx = num2str(beta(i));
    figure()
    subplot(3,1,1)
    plot(1:k+1, x1, '*', 1:k+1, zeros(1, k+1), '-')
    title(strcat('Estimate of x_1 for beta=', idx))

    subplot(3,1,2)
    plot(1:k+1, x2, '*', 1:k+1, zeros(1, k+1), '-')
    title(strcat('Estimate of x_2 for beta=', idx))

    subplot(3,1,3)
    plot(1:k+1, x3, '*', 1:k+1, zeros(1, k+1), '-')
    title(strcat('Estimate of x_3 for beta=', idx))
    hold off
    
    idx = num2str(beta(i));
    figure()
    subplot(3,1,1)
    plot(1:length(fx1), fx1, '*', 1:length(fx1), zeros(1, length(fx1)), '-')
    title(strcat('Estimate of x_1 with the objective function for beta=', idx))

    subplot(3,1,2)
    plot(1:length(fx2), fx2, '*', 1:length(fx2), zeros(1, length(fx2)), '-')
    title(strcat('Estimate of x_2 with the objective function for beta=', idx))

    subplot(3,1,3)
    plot(1:length(fx3), fx3, '*', 1:length(fx3), zeros(1, length(fx3)), '-')
    title(strcat('Estimate of x_3 with the objective function for beta=', idx))
    hold off
end
