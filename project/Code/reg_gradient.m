function [grad] = reg_gradient(num_sensors, Pairwise_Sensor_Distance, Sensor_Anchor_Distance, anchors, x)
%
% Returns the gradient of the regression function
% Inputs:
    % num_sensors: number of sensors to estimate, 
    % anchors: the known anchors
    % Pairwise_Sensor_Distance: sensor1_id, sensor2_id, distance
    % Sensor_Anchor_Distance: sensor_id, anchor_id, distance
    % x: the guess for the sensors!
% return:
	% rad: gradient vector
    
d = size(anchors,1);
grad = zeros(d, num_sensors);


for i=1:size(Pairwise_Sensor_Distance,1)
    sensor1_idx = Pairwise_Sensor_Distance(i,1);
    sensor2_idx = Pairwise_Sensor_Distance(i,2);
    dist = Pairwise_Sensor_Distance(i,3);
    
    x1 = x(:, sensor1_idx);
    x2 = x(:,sensor2_idx);
    diff = x1 - x2;
    
    grad(:, sensor1_idx) = grad(:, sensor1_idx) + 4*(norm(diff)^2 - dist^2)*(diff);
    grad(:, sensor2_idx) = grad(:, sensor2_idx) - 4*(norm(diff)^2 - dist^2)*(diff);
    
end



for i=1:size(Sensor_Anchor_Distance,1)
    sensor_idx = Sensor_Anchor_Distance(i,1);
    anchor_idx = Sensor_Anchor_Distance(i,2);
    dist = Sensor_Anchor_Distance(i,3);
    
    x1 = x(:,sensor_idx);
    anchor = anchors(:,anchor_idx);
    diff = anchor - x1;
    
    grad(:, sensor_idx) = grad(:, sensor_idx) - 4*(norm(diff)^2 - dist^2)*(diff);
    
    
end


end





