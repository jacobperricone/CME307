%% Solves the convex problem with the SDP relaxation.
function [Z] = SDP(num_sensors, Pairwise_Sensor_Distance, Sensor_Anchor_Distance, anchors)

% Number of sensors and anchors.
num_s = length(Pairwise_Sensor_Distance(:, 1));
num_a = length(Sensor_Anchor_Distance(:, 1));

d = length(anchors(:, 1));
n = num_sensors;

cvx_begin sdp quiet
            %% Discuss what do to with the size here.
            variable Z(n+d,n+d) symmetric
            minimize(0)
            subject to
            
            Z(1:d, 1:d) == eye(d, d);
            
            for i=1:num_s
                ei = zeros(1, num_sensors);
                ei(1, Pairwise_Sensor_Distance(i, 1):end)=1;
                ej = zeros(1, num_sensors);
                ej(1, Pairwise_Sensor_Distance(i, 1):end)=1;
                sum(dot([zeros(1, d)';(ei-ej)']*[zeros(1, d)';(ei-ej)']', Z)) == Pairwise_Sensor_Distance(i, 3)^2;
            end
            
            
            for i = 1:num_a
                ej = zeros(1, num_sensors);
                ej(1, Sensor_Anchor_Distance(i, 1):end)=1;
                sum(dot([anchors(:, Sensor_Anchor_Distance(i, 2)); -ej'] * [anchors(:, Sensor_Anchor_Distance(i, 2)); -ej']',  Z))...
                    == Sensor_Anchor_Distance(i, 3)^2;
            end

            Z >= 0;           
cvx_end
        
end 
