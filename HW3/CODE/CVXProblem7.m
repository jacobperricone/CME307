%
% Returns: estimate of the sensor location

function [x] = CVXProblem7( A, b, mu )
cvx_expert true
cvx_begin quiet
    variable x(3);
    opt = .5 * dot(A * x - b, A * x - b) - mu * sum(log(x));

    minimize(opt);
    subject to
                x > 0;
 
cvx_end;
end
