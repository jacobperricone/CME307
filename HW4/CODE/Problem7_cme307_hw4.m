clear all
close all

%% Problem 7
% Part (a) and (b)

A = [1, 1, 1; 1, 1, 2; 1, 2, 2];
beta = [.1, 1, 10];

for i=1:3
    % Without the objective function.
    nA(i,: ) = ADMM(0, beta(i), 199, A);
    
    % With the objective function.
    nAo(i, :) = ADMM(1, beta(i), 199, A);
    
end
figure()
plot(1:200, nA(1, :))
hold on
plot(1:200, nA(2, :))
hold on
plot(1:200, nA(3, :))
title('||Ax|| for the ADMM')
legend('beta=.1', 'beta=1', 'beta=10')
hold off

figure()
plot(1:200, nAo(1, :))
hold on
plot(1:200, nAo(2, :))
hold on
plot(1:200, nAo(3, :))
title('||Ax|| for theADMM with Objective')
legend('beta=.1', 'beta=1', 'beta=10')
hold off



%% Random Permutation
% Part (c)

for i=1:3
    
    % Without the objective function and with permuation.
    nAp(i, :) = PermADMM(0, beta(i), 199, A);
    
    % With objective function and permuation.
    nApo(i, :) = PermADMM(1, beta(i), 199, A);
  
end
figure()
plot(1:200, nAp(1, :))
hold on
plot(1:200, nAp(2, :))
hold on
plot(1:200, nAp(3, :))
title('||Ax|| for the permuated ADMM with Objective')
legend('beta=.1', 'beta=1', 'beta=10')
xlabel('Iterations')
ylabel('||Ax||')

figure()
plot(1:200, nApo(1, :))
hold on
plot(1:200, nApo(2, :))
hold on
plot(1:200, nApo(3, :))
title('||Ax|| for the permuated ADMM with Objective')
legend('beta=.1', 'beta=1', 'beta=10')
xlabel('Iterations')
ylabel('||Ax||')
hold off