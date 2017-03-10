
function [fval] = reg_fval(Pairwise_Sensor_Distance, Sensor_Anchor_Distance, anchors, x)

fval = 0;
for i=1:size(Pairwise_Sensor_Distance,1)
    sensor1_idx = Pairwise_Sensor_Distance(i,1);
    sensor2_idx = Pairwise_Sensor_Distance(i,2);
    dist = Pairwise_Sensor_Distance(i,3);
    
    x1 = x(:, sensor1_idx);
    x2 = x(:,sensor2_idx);
    diff = x1 - x2;
    
    fval = fval + (norm(diff)^2 - dist^2)^2;
    
end



for i=1:size(Sensor_Anchor_Distance,1)
    sensor_idx = Sensor_Anchor_Distance(i,1);
    anchor_idx = Sensor_Anchor_Distance(i,2);
    dist = Sensor_Anchor_Distance(i,3);
    
    x1 = x(:,sensor_idx);
    anchor = anchors(:,anchor_idx);
    diff = anchor - x1;
    
    fval = fval + ( norm(diff)^2 - dist^2)^2;
    
    
end


end




