% question3_client.m
% 
% Dependencies: Anchors.m, VerticleBarPlot.m, and Mesh3D.m
%
% Input: None
%
% This program generates the true and estiamted locations by calling the 
% Anchors and the TDOA functions repsectively. Also generates plots by
% executing the VerticleBarPlot and the Mesh3D.
%
% Output: Plots with locations and errors of estimations.

%%
clear all
close all

tic
%% generates anchors
% Defines dimension of the sensors and anchors.
dim = 1;
anchors = Anchors(dim);
% Outputs a matrix with the true and estimated locations of the sensor.

x_true = 2*rand(1,dim)-.5;

% Fixed sensor, which location is known.
x_tau = 2*rand(1,dim)-.5;
    

% Estimates the location of x_true
estimated_sensors_TDOA = TDOA(anchors, x_true, x_tau);
error_sensors_TDOA = norm(x_true'-estimated_sensors_TDOA);


% verticle bar plots
VerticleBarPlot(estimated_sensors_TDOA, x_true', error_sensors_TDOA, anchors, {' TDOA '})


%%

% 3D mesh
num = 5;
Mesh3D_Plot2D(dim, anchors, num, x_tau)
    
toc