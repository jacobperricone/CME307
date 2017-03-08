function [x_samples] = generate_x_points(npoints)
        x1 = linspace(0.0001, .9999, npoints);
        x2 = rand(size(x1)).*( 1 - x1);
        x3 = 1 - x1 - x2;
        
        x_samples = vertcat(x1,x2,x3);
        

end
