function [ val ] = grad( x, a, d )
val = -4*sum((repmat(dot(a-x, a-x) - d.^2, 2, 1) .* (a-x))')';
end
