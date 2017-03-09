%
% Dependencies: gradx1 and gradx2
%
% Computes the gradient descent for two sensors.


function [x] = SDMwFtwoXs(alpha, a, d, x1_0, x2_0)

x1(:, 1) = x1_0;
x2(:, 1) = x2_0;

for i = 1: 15000 -1 
    x1(:, i+1) = x1(:, i) - alpha * gradx1(a, d, x1(:, i), x2(:, i));
    x2(:, i+1) = x2(:, i) - alpha * gradx2(a, d, x1(:, i+1), x2(:, i));
end

x = [x1(:, end), x2(:, end)];
end
