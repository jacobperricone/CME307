function [  ] = Mesh3D_Plot2D( dim, anchors, num, x_tau )
%% Solves the convex problem with the SDP relaxation.
% Mesh3D.m
%
% Dependencies: TDOA.m
%
% Inputs: dimnsion of the sensor, anchors for the 1D-3D plots, number of
% iterations that generates num by num mesh of error and x_tau, which is
% the fixed sensor for refrence in the TODA function.
%
% This function generates error plots for the 1D-3D TDOA relaxation
% problem. For each iteration it runs the TDOA that estimates the location 
% of a single sensor. 
%
% Returns: Generates a 3D mesh plot (dim = 2,3) of the estimation error, 
% and 2D verticle bar plot of the error when dim = 1. 

%%
    if dim == 1
        x_true = linspace(-2.5, 2.5, num^2)';
    
    elseif dim == 2 
        x_true = linspace(-2.5, 2.5, num)'.*ones(1, dim);
    
    else
        x_true = linspace(-2.5, 2.5, num)'.*[ones(1, dim-1), 0];
    end


    
    %% Generates error mesh or 2D error plots if dim=1.
    if dim == 1

        for i = 1:num^2

            estimated_sensors_TDOA = TDOA(anchors, x_true(i, 1), x_tau);
            error_sensors_TDOA(i) = norm(x_true(i, 1)'- estimated_sensors_TDOA);        
            
        end

        %% Plots of Error
        figure()
        % Convehull
        plot(anchors, [0, 0], 'Color', 'yellow', 'LineWidth', 12) 
        alpha(.4)
        hold on
        plot(x_true, zeros(1, length(x_true)), 'm.', x_true, ...
            error_sensors_TDOA, 'b.-' ,'MarkerSize', 25) 
        hold on
        
        for i = 1: length(error_sensors_TDOA)
            plot([x_true(i), x_true(i)], [0, error_sensors_TDOA(i)], 'k')
            hold on
        end

        legend('Convexhull', 'Error', 'True Location')
        
        xlabel('x')
        ylabel('Error')
        axis([min(x_true)-.5, max(x_true)+.5, 0, max(error_sensors_TDOA)+1])
        [mx,k] = min(error_sensors_TDOA(:));
        [ix,jx] = ind2sub(size(error_sensors_TDOA),k);
        dim = [.10 .595 .3 .3];
        str = strcat('Minimum Error of (', num2str(error_sensors_TDOA(ix,jx)),')', ...
            ' at X = ', num2str(x_true(find(error_sensors_TDOA == min(error_sensors_TDOA)))));
       
          
        
    else    
   
        [X_1, Y_1] = meshgrid(x_true(:, 1), x_true(:, 2));

        for i = 1:num
            for j = 1:num
        
                if dim == 2
                    estimated_sensors_TDOA = TDOA(anchors, [x_true(i, 1), ...
                        x_true(j, 2)], x_tau);
                    error_sensors_TDOA(i, j) = norm([x_true(i, 1), x_true(j, 2)]' ...
                        - estimated_sensors_TDOA);
            
                else 
                    estimated_sensors_TDOA = TDOA(anchors, [x_true(i, 1), ...
                        x_true(j, 2), 0], x_tau);
                    error_sensors_TDOA(i, j) = norm([x_true(i, 1), x_true(j, 2), 0]' ...
                        - estimated_sensors_TDOA);
            
                end
                 
            end
        end

        %% Plots of Error
        figure()
        if dim == 2
            % creates the convexhull around (x,y)
            DT = delaunayTriangulation(anchors');
            % indecies of vertices of the covexhull
            k = convexHull(DT);
        
            plot3(DT.Points(k,1), DT.Points(k,2), 0*DT.Points(k,1),'r')
            legend('Convexhull', 'Location', 'northeast')
            grid on; hold on
            
    
        end
    
        surf(X_1, Y_1, error_sensors_TDOA)
        colormap hsv
        alpha(.4)
        colorbar
        view(-30,30); camlight; axis image

        xlabel('X_1')
        ylabel('X_2')
        zlabel('Error')
        [mx,k] = min(error_sensors_TDOA(:));
        [ix,jx] = ind2sub(size(error_sensors_TDOA),k);
        dim = [.10 .595 .3 .3];
        str = strcat('Minimum Error of (', num2str(error_sensors_TDOA(ix,jx)),')', ...
            ' at X_2: (', num2str(X_1(ix,jx)), ',', num2str(Y_1(ix,jx)),')');
        

    end 
    title('Magntitude of Error for the TDOA', 'FontSize', 10)
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'FontSize',10);

end

