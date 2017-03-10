function [estimated_sensors] = Regression_SDM_Q1(f,df, x_initial, MAX_ITER, TOL, ALPHA, num_sensors,  Pairwise_Sensor_Distance, Sensor_Anchor_Distance, anchors, debug)

    %% find number of unique sensors involved in the optimization
    
    x_prev = x_initial;
    
    iter = 1;
    fvals = [];
    gvals = [];
    
    psd = Pairwise_Sensor_Distance;
    sad = Sensor_Anchor_Distance;
    fvals(iter) = f(psd, sad, anchors, x_prev);
    gvals(iter) = norm(df(num_sensors,psd,sad,anchors,x_prev));
    
    
    
    
    while iter < MAX_ITER
        if abs(gvals(iter)) < TOL
            disp(sprintf('----GRADIENT NORM IS BELOW TOLERANCE CONVERGENCE OF FUNCTION AFTER %d ITERATIONS-----', iter))
            disp(sprintf('-----------------------FINAL Iteration: %d--------------------------------', iter));
            disp('  f(x)    delta_f(x)    g(x)     delta_g(x)   NORMGRAD')
            disp([  fvals(iter), fvals(iter) - fvals(iter-1), gvals(iter), gvals(iter) - gvals(iter) - 1 ])
            break;
        end

            
        iter = iter + 1;
        x_new = x_prev - ALPHA*df(num_sensors, psd, sad, anchors, x_prev);
        
        fvals(iter) = f(psd, sad, anchors, x_new);
        gvals(iter) = norm(df(num_sensors,psd,sad,anchors,x_new));
        
        
        if debug
            disp(sprintf('-----------------------Iteration: %d--------------------------------', iter));
            disp('  f(x)    delta_f(x)    g(x)     delta_g(x)   NORMGRAD')
            disp([  fvals(iter), fvals(iter) - fvals(iter-1), gvals(iter), gvals(iter) - gvals(iter) ])
        end
        
        x_prev = x_new;
    end
        
        
    estimated_sensors = x_prev;
    
    
    
    
    

end