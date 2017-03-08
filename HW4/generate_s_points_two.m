function [spoints] = generate_s_points_two(npoints)
    y = - rand(1,npoints)
    spoints = vertcat(1 - y, -y,1- y)
    
end
