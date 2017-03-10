function [estimated_sensors, fvals, gvals, iter] = Regression_SDM_Q1(f,df, x_initial, MAX_ITER, TOL, ALPHA, num_sensors,  Pairwise_Sensor_Distance, Sensor_Anchor_Distance, anchors, debug)
%
% Description: Implements the Steepest Descent Method on an initial point
% x_initial, generates plots of function value and gradient value over time
% Inputs
    % f: function to optimize over-- takes arguments in the form
    %    f(Pairwise_Sensor_Distance, Sensor_Anchor_Distance, anchors, x)
    %
    % df: gradient function that takes arguments in the form
    %     df(num_sensors, Pairwise_Sensor_Distance, Sensor_Anchor_Distance, anchors, x)
    %     where x is the point at which to calculate the gradient 
    %
    % num_sensors: number of sensors to estimate, 
    % anchors: the known anchors
    % Pairwise_Sensor_Distance: sensor1_id, sensor2_id, distance
    % Sensor_Anchor_Distance: sensor_id, anchor_id, distance
    % MAX_ITER, TOL, ALPHA, debug all kind of obvious
% Returns: 
%   estimated_sensors: estimated sensors 
    
    x_prev = x_initial;
    
    iter = 1;
    fvals = [];
    gvals = [];
    
    psd = Pairwise_Sensor_Distance;
    sad = Sensor_Anchor_Distance;
    fvals(iter) = f(psd, sad, anchors, x_prev);
    df(num_sensors,psd,sad,anchors,x_prev)
    gvals(iter) = norm(df(num_sensors,psd,sad,anchors,x_prev));
    
    
    
    
    while iter < MAX_ITER
        if abs(gvals(iter)) < TOL
            disp(sprintf('----GRADIENT NORM IS BELOW TOLERANCE CONVERGENCE OF FUNCTION AFTER %d ITERATIONS-----', iter))
            disp(sprintf('-----------------------FINAL Iteration: %d--------------------------------', iter));
            disp('  f(x)    delta_f(x)    g(x)     delta_g(x)   NORMGRAD')
            disp([  fvals(iter), fvals(iter) - fvals(iter-1), gvals(iter), gvals(iter) - gvals(iter)])
            break;
        end

            
        iter = iter + 1;
        x_new = x_prev - ALPHA*df(num_sensors, psd, sad, anchors, x_prev);
        
        fvals(iter) = f(psd, sad, anchors, x_new);
        gvals(iter) = norm(df(num_sensors,psd,sad,anchors,x_new));
        
        
        if debug
            disp(sprintf('-----------------------Iteration: %d--------------------------------', iter));
            disp('  f(x)       delta_f(x)       g(x)        delta_g(x)      NORMGRAD')
            disp([  fvals(iter), fvals(iter) - fvals(iter-1), gvals(iter), gvals(iter) - gvals(iter) ])
        end
        
        x_prev = x_new;
    end
        
        
    estimated_sensors = x_prev;
    
    
    
    
s = 'Regression';
fig = figure();
subplot(2,1,1)
plot(1:iter, gvals, 'LineWidth',2); grid on;
title(strcat({' Norm Gradient Function SDM Method: '}, s));
xlabel('Iteration'); ylabel('G(x)');

subplot(2, 1, 2)
plot(1:iter, fvals, 'LineWidth',2); grid on;
title(strcat({'Objective Function SDM Method: '}, s));

xlabel('Iteration'); ylabel('F(x)');

save_string = [pwd strcat('/FIGURES/Question1/SDM', '_', s,'.png')];
saveas(fig, save_string)


end
