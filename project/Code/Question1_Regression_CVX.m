function [estimated_sensors] = SDP(num_sensors, Pairwise_Sensor_Distance, Sensor_Anchor_Distance, anchors)

% Number of sensors and anchors iterations.
num_s = length(Pairwise_Sensor_Distance(:, 1));
num_a = length(Sensor_Anchor_Distance(:, 1));

% Number of sensors and anchors
d = length(anchors(:, 1));
n = num_sensors;


% Number of sensors and anchors
d = length(anchors(:, 1));
n = num_sensors;


cvx_begin quiet
            %% Discuss what do to with the size here.
            
            variable x(d, n)
            variable obj
            
            for i=1:num_s
                obj = obj + (norm(x(:,Pairwise_Sensor_Distance(i,1)) - x(:,Pairwise_Sensor_Distance(i,2)), 'fro')^2 - Pairwise_Sensor_Distance(i,3)^2)^2
            end
            
            for i=1:num_a
                obj = obj + (norm(Sensor_Anchor_Distance(i,2) - x(:,Pairwise_Sensor_Distance(i,2)), 'fro')^2 - Sensor_Anchor_Distance(i,3)^2)^2
            end
            
            minimize(obj)
            
cvx_end

estimated_sensors = x



end
