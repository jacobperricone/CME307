function [Z] = SDP_Noise(num_sensors, Pairwise_Sensor_Distance, Sensor_Anchor_Distance, anchors)

% Number of sensors and anchors iterations.
num_s = length(Pairwise_Sensor_Distance(:, 1));
num_a = length(Sensor_Anchor_Distance(:, 1));


% Number of sensors and anchors
d = length(anchors(:, 1));
n = num_sensors;


cvx_begin sdp quiet
            %% Discuss what do to with the size here.
            variable Z(n+d, n+d) symmetric;
            variable u_ij(num_s);
            variable v_ij(num_s);
            variable u_ik(num_a);
            variable v_ik(num_a);
            z = 0;
            for i=1:num_s
                z = z + u_ij(i)+ v_ij(num_s);
            end
            
            for i=1:num_a
                z = z +u_ik(i) +v_ik(i);
            end
            
            
            minimize(z);
            
            subject to
            
           Z(1:d, 1:d) == eye(d, d);
            
            for i=1:num_s
                ei = zeros(1, num_sensors);
                ei(1, Pairwise_Sensor_Distance(i, 1)) = 1;
                
                ej = zeros(1, num_sensors);
                ej(1, Pairwise_Sensor_Distance(i, 2)) = 1;
                
                sum(dot([zeros(1, d)';(ei-ej)']*[zeros(1, d)';(ei-ej)']', Z)) - u_ij(i) + v_ij(i) == Pairwise_Sensor_Distance(i, 3)^2;
            end
            
            
            for j = 1:num_a
                ej = zeros(1, num_sensors);
                ej(1, Sensor_Anchor_Distance(j, 1)) = 1;
                
                sum(dot([anchors(:, int64(Sensor_Anchor_Distance(j, 2))); -ej'] * [anchors(:, int64(Sensor_Anchor_Distance(j, 2))); -ej']',  Z))...
                    - u_ik(j) + v_ik(j)== Sensor_Anchor_Distance(j, 3)^2;
            end

            Z >= 0; 
            u_ij >= 0;
            v_ij >= 0;
            u_ik >= 0;
            v_ik >= 0;
            
            
cvx_end


end