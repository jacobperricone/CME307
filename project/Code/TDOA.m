function [output] = TDOA(anchors, x_true, x_tau)
%% Solves the convex problem with the SDP relaxation.
% TDOA.m
%
% Inputs: anchors.
%
% This funciton estimates locations of n sensors in 2 or 3D plane by 
% using CVX.
%
% Returns: dim by 2 matrix with first column corresponding to the true
% sensor location while the second is the estimte given by the cvx.

    d = size(anchors, 1);
    % Locations of the anchors. Each enchor corresponds to the column of a.
    anchors = Anchors(d);
    num_anchors = size(anchors, 2);
    
    anchor_dist = pdist2(x_true, anchors');
    d_tau = pdist2(x_true, x_tau);

    %%
    cvx_expert true
    cvx_begin sdp quiet

        variable t
        variable ys 
        variable ds
        variable Y(d+1, d+1) symmetric
            
        minimize(t)
        subject to
            
            Y(1:d, 1:d) == eye(d);
                          
            for i=1:num_anchors
            
                X = zeros(d+1, d+1);    
                X(1:d, 1:d) = kron(anchors(:, i), anchors(:, i)');
                X(end, 1:d) = -anchors(:, i)';
                X(1:d, end) = -anchors(:, i);
                X(end, end) = 1;
                
                % Given that cvx doesnt distinguish between < and <= I
                % replace < with <= so it doesn't throw warnings.
                
                trace(Y*X) - trace([1, d_tau; d_tau, ds]*[anchor_dist(i)^2, ...
                    anchor_dist(i); anchor_dist(i), 1]) >= - t
                trace(Y*X) - trace([1, d_tau; d_tau, ds]*[anchor_dist(i)^2, ...
                    anchor_dist(i); anchor_dist(i), 1]) <= t    
            end
             
            
            ds >= d_tau^2;
            ds == dot(x_tau, x_tau) - 2*dot(x_tau, Y(end, 1:d)') + ys;
            X >= 0;
            Y >= 0
            Y(end, end) == ys;
            ys >= dot(Y(end, 1:d), Y(end, 1:d));
            [1, d_tau; d_tau, ds] >= 0
                     
    cvx_end

    %%
    output = Y(end, 1:d)';

end