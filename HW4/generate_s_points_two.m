function [spoints] = generate_s_points_two(npoints)
    y =  rand(1,npoints);
    spoints = vertcat(1 + y, y, y);
    
end

