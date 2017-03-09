clear all
close all


%% generates anchors
random_x = [3*rand(1,6)-1.5]';
random_y = [3*rand(1,6)-1.5]';
% creates the convexhull around (x,y)
DT = delaunayTriangulation(random_x, random_y);
% indecies of vertices of the covexhull
k = convexHull(DT);
anchors = [(3*DT.Points(k(1:end-1),1)-1)'; (3*DT.Points(k(1:end-1),2)-1)'];
size_x = size(anchors,1);
num_anchors = size(anchors,2);
% find number of sensors 
num_sensors = 40;
% generate them
[sensors, dx, da] = generate_sensor(anchors,num_sensors,1, DT, k);

% Go from index to x,y coordinat
id_xy = @(i,rows, cols) ([ (mod(i,rows) == 0)*rows + mod(i,rows)  , ceil(i/rows)]);
% go from x,y to index
id = @(row,col,rows, columns) ((col-1)*rows + row);

%% get the indices of the upper quadrant of the d_x matrix

upper_quadrant = zeros(1, round(num_sensors^2/2) - num_sensors);
k = 1;
for i=1:num_sensors
    for j=i+1:num_sensors
        upper_quadrant(k) = id(i,j, num_sensors, num_sensors);
        k = k + 1;
    end
end

% Number of sensors for which we have pairwise distance information
num_Nx = round(numel(upper_quadrant));
% Index of sensors for which we have anchor information
Nx = upper_quadrant(randperm(numel(upper_quadrant), num_Nx))';
% find x,y distance
X_x = id_xy(Nx, num_sensors, num_sensors);
% create matrix

tmpX = zeros(size(X_x,1),1);
for i=1:size(X_x,1)
    tmpX(i) = dx(X_x(i,1), X_x(i,2));
end


Pairwise_Sensor_Distance = horzcat(X_x, dx(Nx));
%%
% -2*num_sensors works only for num_sensors>5, otherwise it needs to be
% smaller.
num_Na = round(numel(da))-2.5*num_sensors;
Na = randperm(numel(da), num_Na)';
X_a = id_xy(Na, num_sensors, num_anchors);

tmpA = zeros(size(X_a,1),1);
for i=1:size(X_a,1)
    tmpA(i) = da(X_a(i,1), X_a(i,2));
end

Sensor_Anchor_Distance = horzcat(X_a, da(Na));

%% solves the SDP Relaxation Problem.
Z = SDP(num_sensors, Pairwise_Sensor_Distance, Sensor_Anchor_Distance, anchors);
estimated_sensors = Z(length(anchors(:,1))+1:end, 1:size_x)';

error_sensors = ((sensors(1, :) - estimated_sensors(1, :)).^2 + ...
    (sensors(2, :) - estimated_sensors(2, :)).^2).^.5;


%% 3d plot
ThreeDVerticleBarPlot(estimated_sensors, sensors, error_sensors)