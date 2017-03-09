%%
% ThreeDVerticleBarPlot
% Input: vectors of true sensor locations, estiamted sensor locations, and 
% the magnitude of the error that comes from the estimation using CVX. All 
% sensor locations are defined in 2D. 
% Retuns: 3D bar plot. Adds estiamted sensor locations to the 2D plot 
% initiatated in generate_sensor.m

function [  ] = ThreeDVerticleBarPlot(estimated_sensors, sensors, error_sensors)

    % Add estiamted locations to the 2D plot.
    hold on
    scatter(estimated_sensors(1,:), estimated_sensors(2,:), 'm.')

    % 3D plot
    figure()
    x = sensors(1, :);
    y = sensors(2, :);
    z = zeros(1, length(y));
    h = plot3(sensors(1, :), sensors(2,:), zeros(1, length(sensors(1,:))), '.r');
    set(h, 'MarkerSize', 15);
    grid on;  hold on

    
    for i=1:length(sensors(1, :))
        x_verticle = [x(i); x(i)];
        y_verticle = [y(i); y(i)];
        z_verticle = [z(i) + error_sensors(i), z(i)];
        
        % Draws vertical errorbar for each point.
        h=plot3(x_verticle, y_verticle, z_verticle, '-k');
        set(h, 'LineWidth', 1.5, 'Color', 'b');
    end
    
    title('Magnitude of the Error for Each Estiamtion of the Sensor Location.')
    xlabel('x_1')
    ylabel('y_1')
    zlabel('Error')

end

