clear all
close all

%% Problem 7
% Part (a) and (b)

A = [1, 1, 1; 1, 1, 2; 1, 2, 2];
beta = [.1, 1, 10];

for i=1:3
    % Without the objective function.
    ADMM(0, beta(i), 20, A)
    
    % With the objective function.
    ADMM(1, beta(i), 20, A)
end




%% Random Permutation
% Part (c)

for i=1:3
    
    % Without the objective function and with permuation.
    PermADMM(0, beta(i), 20, A) 
    
    % With objective function and permuation.
    PermADMM(1, beta(i), 20, A) 
  
end
