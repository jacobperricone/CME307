function [sensors, d_x, d_a] = generate_sensor(anchors, num_sensors, debug)
%
% generate_sensors
% inputs: anchor points that specify the convex hull of the region (d x N),
    % where d is the number of dimesnions, and N is the number of anchors
    % num_sensors: number of sensors to generate 
    % returns: num_sensors randomly generate sensors
    % d_x: distance between all sensors
    % d_a: distance between all sensors and all anchors
% returns:
    % sensors: randomly generate sensors
    % d_x: distance between the sensors
    % d_a: distance between anchors and sensors
    
    num_dimensions = size(anchors,1);
    min_values = min(anchors,[],2);
    max_values = max(anchors,[],2);
    
    sensors = repmat(min_values  - mean(max_values - min_values,2)/4, 1,num_sensors,1)...
        + rand(num_dimensions,num_sensors).*(repmat(max_values + mean(max_values - min_values,2)/2 - min_values, 1, num_sensors,1));
     
    
    
    
    if debug
        
        if num_dimensions == 1
            plot(sensors, zeros(1, length(sensors)), 'ko') 
            
        elseif num_dimensions == 2
            scatter(sensors(1,:), sensors(2,:), 'ko')
            
        elseif num_dimensions == 3
            plot3(sensors(1,:), sensors(2,:), sensors(3,:), 'ko')
            
        end
        
    end
   
    
    
    d_x = pdist2(sensors',sensors');
    d_a = pdist2(sensors',anchors');
    
    

end
