%%
clear all
close all
a = [1, -1, 0; 0, 0, 2];
alpha = .001;

for k=1:2
    [X_1, Y_1] = meshgrid(linspace(-1.5, 1.5,5)',linspace(-.5,2.5,5)');
    norm_SPD = zeros(size(X_1));
    norm_SPD_Rand= zeros(size(X_1));
%%
    noise = rand(1) - 1;
    A= [1, 0, 1; 0, 1,1; 0,0,0] ;
    for j=1:size(X_1, 1)
        for l=1:size(Y_1,1)
        
            d = pdist2([X_1(j,l),Y_1(j,l)], a');
        
            cvx_begin sdp quiet
                variable Z(3,3) symmetric
                minimize(0);
                subject to

                sum(dot(A(:,1)*A(:,1)', Z)) == 1;
                sum(dot(A(:,2)*A(:,2)', Z)) == 1;
                sum(dot(A(:,3)*A(:,3)', Z)) == 2;

                for i = 1:3
                    sum(dot([a(:, i); -1] * [a(:, i); -1]',  Z)) == d(i)^2;
                end

                Z >= 0;

            cvx_end
        
            if k == 2
                d = d + noise;
                idx = strcat(' with d Noise of ', num2str(noise));
            elseif k == 1
                idx = ' without d Noise';
            end
        
            x = SDMwF(alpha, a, d, Z(end,1:2)');
            norm_SPD(j,l) = norm(x - [X_1(j,l), Y_1(j,l)]);
        
        
            xrand = [3*rand(1)-1.5, 3*rand(1)-.5];
            x = SDMwF(alpha, a, d, xrand');
            norm_SPD_Rand(j,l) = norm(x - [X_1(j,l), Y_1(j,l)]);
        
    
        end
    
    end

%%
    figure()
    mesh(X_1, Y_1, norm_SPD)
    title(strcat('Magntitude of Error with CVX Computed Starting Point x_0', idx))
    xlabel('X Data')
    ylabel('Y Data')
    zlabel('Error')

    figure()
    mesh(X_1, Y_1, norm_SPD_Rand)
    title(strcat('Magntitude of Error with Random Starting Point x_0', idx))
    xlabel('X Data')
    ylabel('Y Data')
    zlabel('Error')
    hold off
end
