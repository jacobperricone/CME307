%%
% Client Script For Problem two
%
%%
clear all
close all

tic
%% generates anchors, and distance information
dim = 1;
anchors = Anchors(dim);



size_x = size(anchors, 1);
num_anchors = size(anchors,2);
% find number of sensors 
num_sensors = 40;
% generate them
[sensors, dx, da] = generate_sensor(anchors,num_sensors,1);

% Go from index to x,y coordinat
id_xy = @(i,rows, cols) ([ (mod(i,rows) == 0)*rows + mod(i,rows)  , ceil(i/rows)]);
% go from x,y to index
id = @(row,col,rows, columns) ((col-1)*rows + row);

% get the indices of the upper quadrant of the d_x matrix

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
%
if dim == 1
    num_Na = round(numel(da))-num_sensors;
else
    num_Na = round(numel(da))-2.5*num_sensors;
end

Na = randperm(numel(da), num_Na)';
X_a = id_xy(Na, num_sensors, num_anchors);

tmpA = zeros(size(X_a,1),1);
for i=1:size(X_a,1)
    tmpA(i) = da(X_a(i,1), X_a(i,2));
end

Sensor_Anchor_Distance = horzcat(X_a, da(Na));
%% 
Z = SDP_Noise(num_sensors, Pairwise_Sensor_Distance, Sensor_Anchor_Distance, anchors);
estimated_sensors_SDP = Z(length(anchors(:,1))+1:end, 1:size_x)';



if dim == 1
    error_sensors_SDP = ((sensors'-estimated_sensors_SDP').^2).^.5;
else 
    error_sensors_SDP = sum(((sensors-estimated_sensors_SDP).^2)).^.5;
end

% verticle bar plots
VerticleBarPlot(estimated_sensors_SDP, sensors, error_sensors_SDP, anchors, 'SDP No Descent')
%% Add SDM

df = @reg_gradient 
f = @reg_fval 
% Set the convergence and iteration limits
MAX_ITER = 1000; TOL = .0001; ALPHA = .0001; debug = 0; 
% run the descent methods

%%
[estimated_sensors_SDM, fvals_SDM, gvals_SDM,iter_SDM] = Regression_SDM_Q1(f,df, estimated_sensors_SDP, MAX_ITER, TOL, ALPHA, num_sensors,  Pairwise_Sensor_Distance, Sensor_Anchor_Distance, anchors, debug);
%%
[estimated_sensors_ASDM, fvals_ASDM, gvals_ASDM, iter_ASDM]= Regression_ASDM_Q1(f,df, estimated_sensors_SDP, MAX_ITER, TOL, ALPHA, num_sensors,  Pairwise_Sensor_Distance, Sensor_Anchor_Distance, anchors, debug);

if dim == 1
    error_sensors_SDM = ((sensors'-estimated_sensors_SDM').^2).^.5;
    error_sensors_ASDM = ((sensors'-estimated_sensors_ASDM').^2).^.5;
else   
    error_sensors_SDM = sum(((sensors-estimated_sensors_SDM).^2)).^.5;
    error_sensors_ASDM = sum(((sensors-estimated_sensors_ASDM).^2)).^.5;
end
%%
VerticleBarPlot(estimated_sensors_SDM, sensors, error_sensors_SDM,anchors, 'SDP with SDM')
VerticleBarPlot(estimated_sensors_ASDM, sensors, error_sensors_ASDM,anchors,'SDP with ASDM')
%%
figure()
subplot(2,1,1)
loglog(1:iter_SDM, fvals_SDM)
hold on
loglog(1:iter_ASDM, fvals_ASDM)
legend('SDM', 'ASDM')
title('Function Value vs. Iteration ')
xlabel('$\log[\mbox{Iteration}]$','Interpreter', 'LaTex')
ylabel('$\log[f(x)]$', 'Interpreter', 'LaTex')

subplot(2,1,2)
loglog(1:iter_SDM, gvals_SDM)
hold on
loglog(1:iter_ASDM, gvals_ASDM)
legend('SDM', 'ASDM')
title('Norm Gradient Value vs. Iteration ')
xlabel('$\log[\mbox{Iteration}]$','Interpreter', 'LaTex')
ylabel('$\log[g(x)]$', 'Interpreter', 'LaTex')



