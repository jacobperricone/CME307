function [ val ] = grad( x, a, d )
val = sum( repmat(- 4*(diag((a - repmat([1;0], 1, size(a,2)))'*(a - repmat([1;0], 1, size(a,2)))) - b.^2)', 2, 1).*(a - repmat([1;0], 1, size(a,2))),2)
end