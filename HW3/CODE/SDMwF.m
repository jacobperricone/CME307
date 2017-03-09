% SDMwF.m
%
% Dependencies: grad.m
%
% Computes the gradient descent estimate for a single sensor x in 2D
%
% Returns: the estimated location of a sensor using the gradeint descent.


function [x] = SDMwF(alpha, a, d, x0)

xk(:, 1) = x0;
out(1, :) =[xk(:, 1)'];

for i = 1: 15000 - 1 
    xk(:, i+1) = xk(:, i) - alpha * grad(xk(:, i), a, d);
    out(i+1, :) = [xk(:, i+1)'];

end
x = out(end, 1:2);
end
