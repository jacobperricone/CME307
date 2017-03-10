%%
% VerticleBarPlot.m 
% Client: setup_problem_client.m
%
% Input: vectors of true sensor locations, estiamted sensor locations, and 
% the magnitude of the error that comes from the estimation using CVX. All 
% sensor locations are defined in 2D. 
%
% Retuns: 3D bar plot. Adds estiamted sensor locations to the 2D plot 
% initiatated in generate_sensor.m

function [  ] = VerticleBarPlot(estimated_sensors, sensors, error_sensors)

    dim = size(sensors, 1);
    
    if dim == 1
        % Add estiamted locations to the 2D plot.
        hold on
        plot(estimated_sensors, zeros(1, length(estimated_sensors)), 'm.')
        title(strcat('Ture and Estimated Sensor Locations with n=',  ...
            num2str(length(estimated_sensors))))
        xlabel('Interval')
        legend('Convexhull', 'True Locations', 'Estimated Locations')
        
        figure()
        plot(sensors, zeros(1, length(sensors)), '.' ,'MarkerSize', 25) 
        hold on
        for i = 1: length(error_sensors)
            plot([sensors(i), sensors(i)], [0, error_sensors(i)], 'r')
            hold on
        end
        
        title(strcat('Magnitude of the Error for Each Estiamtion of the Sensor Location with n=',  ...
            num2str(length(estimated_sensors))))
        xlabel('x')
        ylabel('Error')
        
        
    else
        % Add estiamted locations to the 2D plot.
        hold on
        if dim == 2
            scatter(estimated_sensors(1,:), estimated_sensors(2,:), 'm.')
            title(strcat('Ture and Estimated Sensor Locations with n=',  ...
                num2str(length(estimated_sensors))))
            %xlabel('x_1')
            %ylabel('x_2')
            %legend('yeah', 'yeah1', 'yeah2')
        else
            plot3(estimated_sensors(1,:), estimated_sensors(2,:), estimated_sensors(3,:), 'm.')
            zlabel('x_3')
            
        end
       
        title(strcat('Ture and Estimated Sensor Locations with n=', ...
            num2str(length(estimated_sensors))))
        xlabel('x_1')
        ylabel('x_2')
        legend('Convexhull', 'True Locations', 'Estimated Locations', 'Location','northeast')

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
        
        
        if dim == 2
            title(strcat('Magnitude of the Error for Each Estiamtion of the Sensor Location with n=', ...
                num2str(length(estimated_sensors))))
        else
            title(strcat('Magnitude of the Error for Each Estiamtion of the Sensor Location with Fixed value of z with n=',...
                num2str(length(estimated_sensors))))
        end
            xlabel('x_1')
            ylabel('y_1')
            zlabel('Error')
      
    end

end
