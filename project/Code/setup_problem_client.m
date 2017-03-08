clear all
close all


%% declare the anchors
anchors = [1, -1, 0; 0, 0, 2];
num_anchors = size(anchors,2);
% find number of sensors 
num_sensors = 10;
% generate them
[sensors, dx, da] = generate_sensor(anchors,num_sensors,0);

% Go from index to x,y coordinat
id_xy = @(i,rows, cols) ([mod(i,rows), floor(i/rows) + 1]);
% go from x,y to index
id = @(row,col,rows, columns) ((col-1)*rows + row);

%% get the indices of the upper quadrant of the d_x matrix
upper_quadrant = zeros(1, num_sensors^2/2 - num_sensors);
k = 1;
for i=1:num_sensors
    for j=i+1:num_sensors
        upper_quadrant(k) = id(i,j, num_sensors, num_sensors);
        k = k + 1;
    end
end

% Number of sensors for which we have pairwise distance information
num_Nx = round(numel(upper_quadrant)/4);
% Index of sensors for which we have anchor information
Nx = upper_quadrant(randperm(numel(upper_quadrant), num_Nx))';
% find x,y distance
X_x = id_xy(Nx, num_sensors, num_sensors);
% create matrix

tmpX = zeros(size(X_x,1));
for i=1:size(X_x,1)
    tmpX(i) = dx(X_x(i,1), X_x(i,2));
end


Pairwise_Sensor_Distance = horzcat(X_x, dx(Nx));


num_Na = round(numel(da)/4);
Na = randperm(numel(da), num_Na)';
X_a = id_xy(Na, num_sensors, num_anchors);

tmpA = zeros(size(X_a,1));
for i=1:size(X_a,1)
    tmpA(i) = da(X_a(i,1), X_a(i,2));
end


Sensor_Anchor_Distance = horzcat(X_a, da(Na));