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

function [  ] = VerticleBarPlot(estimated_sensors, sensors, error_sensors, anchors, name)

    dim = size(sensors, 1);
        
    if dim == 1
        figure()
        % Add estiamted locations to the 2D plot.
        plot(anchors, [0, 0], 'Color', 'cyan', 'LineWidth', 10) 
        axis([anchors(1)-.5, anchors(2)+.5, -.01, .01])
        
        hold on
        plot(sensors, zeros(1, length(sensors)), 'ko')
        hold on
        plot(estimated_sensors, zeros(1, length(estimated_sensors)), 'm.')
        title(strcat({'True and Estimated Sensor Locations Using the ' }, name, {' Method and n='},  ...
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
        
        title(strcat({'||Error|| for Each Sensor Estimate using the '}, name, {' Method and n='},  ...
            num2str(length(estimated_sensors))))
        xlabel('x')
        ylabel('Error')
        
        
    else
        
        % Add estiamted locations to the 2D plot.   
        if dim == 2
            % creates the convexhull around (x,y)
            DT = delaunayTriangulation(anchors');
            % indecies of vertices of the covexhull
            k = convexHull(DT);
        
            figure()
            plot(3*DT.Points(k,1)-1, 3*DT.Points(k,2)-1,'r') 
            hold on
            scatter(sensors(1,:), sensors(2,:), 'ko') 
            hold on
            scatter(estimated_sensors(1,:), estimated_sensors(2,:), 'm.')
           
        else
            % creates the convexhull around (x,y)
            DT = delaunayTriangulation(anchors');
            % indecies of vertices of the covexhull
            [k,v] = convexHull(DT);
        
            figure()
            trisurf(k,3*DT.Points(:,1)-1, 3*DT.Points(:,2)-1, 3*DT.Points(:,3)-1, 'FaceColor', 'cyan')
            alpha(.1)
            hold on 
            plot3(sensors(1,:), sensors(2,:), sensors(3,:), 'ko')
            hold on
            plot3(estimated_sensors(1,:), estimated_sensors(2,:), estimated_sensors(3,:), 'm.')
            zlabel('x_3')
            
        end
       
        title(strcat({'True vs. Estimated Sensor Locations using '}, name, {' Method and n='}, ...
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
            title(strcat({'||Error|| for Each Sensor Estimate using the '}, name, {' Method and n='}, ...
                num2str(length(estimated_sensors))))
        else
            title(strcat({'||Error|| for Each Sensor Estimate using the ' }, name, {' Method and n='},...
                num2str(length(estimated_sensors))))
        end
            xlabel('x_1')
            ylabel('y_1')
            zlabel('Error')
      
    end

end

