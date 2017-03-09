%
% Returns: 2x1 vector corresponding to the gradent of x2.

function [val] = gradx2(a, d, x1, x2)
val = -4*(norm(a(:, 2) - x2)^2 - d(3)^2) * (a(:, 2) - x2) ....
    -4*(norm(a(:, 3) - x2)^2 - d(4)^2) * (a(:, 3) - x2) ...
    -4*(norm(x1 - x2)^2 - d(5)^2)*(x1-x2);
end

