%%
% Problem7_Hw3.m
%
% This program finds the sensor location for a 2d problem by first
% calling the CVXProblem7 function that estimates the location of the
% sensor. It then uses that estimate as an initial point for two sets
% of gradient descent fomulations (general form and the scaled one) in
% the double for-loop.
%
% Retuns: Generates a set of plots


clear all
close all

%% 7(d)
A = [1 0 2; 3 1 2; 1 0 1];
b = [2; 4; 3];
e = [1; 1; 1];
ALPHA = .001;
mu = [0, 10^-4];

for k = 1:2
    xtrue = CVXProblem7(A, b, mu(k));

    f = @(x) .5*norm(A*x - b)^2 - mu(k)*sum(log(x));
    gradL = @(x) A' * (A * x - b)- mu(k) * diag(1 ./ x) * e;
    gradLsc = @(x) x.^2 .* (A' * (A * x - b)) - mu(k) .* x;

    [X_1, Y_1] = meshgrid(linspace(1, 10,30)',linspace(1,10,30)');
    errorL = zeros(size(X_1,1), size(Y_1,1));
    errorLsc = zeros(size(X_1,1), size(Y_1,1));

    for j=1:size(X_1, 1)
        for l=1:size(Y_1,1)

            x(:, 1) = [j; l; 1];
            xsc(:, 1) = [j; l; 1];

            for i = 1:20000
                x(:, i + 1) = x(:, i) - ALPHA .* (gradL(x(:, i)));
                xsc(:, i + 1) = xsc(:, i) - ALPHA .* (gradLsc(xsc(:, i)));
            
            end
            fL(j, l) = norm(f(x(:, end)));
            fLsc(j, l) = norm(f(xsc(:, end)));
            errorL(j, l) = norm(xtrue - x(:, end));
            errorLsc(j, l) = norm(xtrue - xsc(:, end));
        end
    end

    figure()
    subplot(4,1,1)
    mesh(X_1, Y_1, errorL)
    colormap hsv
    alpha(.4)
    colorbar
    view(-30,30); camlight; axis image

    title(strcat('Magntitude of x''s Error for the Standrad SGD with mu=', num2str(mu(k))) , 'FontSize', 10)
    xlabel('Value of x_1')
    ylabel('Value of x_2')
    zlabel('Error')

    
    subplot(4,1,2)
    mesh(X_1, Y_1, fL)
    colormap hsv
    alpha(.4)
    colorbar
    view(-30,30); camlight; axis image

    title('Value of the Objective Function for the Standrad SGD', 'FontSize', 10)
    xlabel('Value of x_1')
    ylabel('Value of x_2')
    zlabel('Cost')
    
    
    subplot(4,1,3)
    mesh(X_1, Y_1, errorLsc)
    colormap hsv
    alpha(.4)
    colorbar
    view(-30,30); camlight; axis image

    title(strcat('Magntitude of x''s Error for the Scaled SGD with mu=', num2str(mu(k))), 'FontSize', 10)
    xlabel('Value of x_1')
    ylabel('Value of x_2')
    zlabel('Error')
    
    
    subplot(4,1,4)
    mesh(X_1, Y_1, fLsc)
    colormap hsv
    alpha(.4)
    colorbar
    view(-30,30); camlight; axis image

    title('Value of the Objective Function for the Scaled SGD', 'FontSize', 10)
    xlabel('Value of x_1')
    ylabel('Value of x_2')
    zlabel('Cost')
    
    hold off
end
