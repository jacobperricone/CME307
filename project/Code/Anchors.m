% Client: setup_problem_client.m
% 
% Input: dimension of the anchor.
%
% The function produces a set of random anchors for a given dim
% as well as generates a plot of the convexhull.
%
% Return: set of anchor locations:
%   dim = 1: two anchors that mark an interval in 1D
%   dim = 2: set of 2D anchors
%   dim = 3: set of 3D anchors




function [ anchors ] = Anchors(dim)

    random_x = [3*rand(1,5)-1.5]';
    random_y = [3*rand(1,5)-1.5]';
    random_z = [3*rand(1,5)-1.5]';

    if dim == 1
        
        anchors = [min(random_x), max(random_x)];

    elseif dim == 2
        
        % creates the convexhull around (x,y)
        DT = delaunayTriangulation(random_x, random_y);
        % indecies of vertices of the covexhull
        k = convexHull(DT);
        anchors = [(3*DT.Points(k(1:end-1),1)-1)'; (3*DT.Points(k(1:end-1),2)-1)'];
        

    elseif dim == 3
        
        % creates the convexhull around (x,y)
        DT = delaunayTriangulation(random_x, random_y, random_z);
        anchors = DT.Points';
           
    end
    
end


