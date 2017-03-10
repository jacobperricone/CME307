% SOCP.m
% Inputs: The total number of sensors num_sensors, a 3 column matrix 
% Pairwise_Sensor_Distance with indecies of the sensor(i), sensor(j) and 
% the distance between them in its columnsm, a 3 column matrix 
% Sensor_Anchor_Distance with columns: sensor indecies, anchor indecies and
% the distance between them, and finally all the anchors.
%
% This funciton estimates locations of n sensors in 1, 2 or 3D plane by 
% using CVX.
%
% Returns: nxd matrix (where x<=d) with all estimates of the sensor 
% locations

%% Solves the convex problem with the SDP relaxation.
function [X] = SOCP(num_sensors, Pairwise_Sensor_Distance, Sensor_Anchor_Distance, anchors)
a = anchors;
% Number of sensors and anchors iterations.
num_s = length(Pairwise_Sensor_Distance(:, 1));
num_a = length(Sensor_Anchor_Distance(:, 1));

% Number of sensors and anchors
d = length(anchors(:, 1));
n = num_sensors;


cvx_begin sdp quiet
            %% Discuss what do to with the size here.
            variable X(d, n)
            minimize(0)
            subject to
            
            
            for i=1:num_s                
                norm(X(:, Pairwise_Sensor_Distance(i, 1))- X(:, Pairwise_Sensor_Distance(i, 2))) <= Pairwise_Sensor_Distance(i, 3);
            end
            
            
            for j = 1:num_a 
                norm(anchors(:, Sensor_Anchor_Distance(j, 2)) - X(:, Sensor_Anchor_Distance(j, 1))) <= Sensor_Anchor_Distance(j, 3);
            end
         
cvx_end
    
end 
