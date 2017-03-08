function [spoints] = generate_s_points(npoints)
    y = - rand(1,npoints);
    spoints = vertcat(1 - y, 1- y, -y);
    
end
