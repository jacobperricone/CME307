%% Inside of the Convex Hull
a = [1, -1, 0; 0, 0, 2];
lambda = rand(3,1);
lambda = lambda / sum(lambda);
xtrue = a * lambda;
d = pdist2(xtrue', a');

cvx_begin quiet
variable x(2);
minimize([0.0]);
subject to
for i = 1:3
    norm(x - a(:, i)) <= d(i);
end
cvx_end
norm(x-xtrue)



%% Outside of the Convex Hull
a = [1, -1, 0; 0, 0, 2];
%lambda = rand(3,1);
%lambda = lambda / sum(lambda);
xtrue = [6; 5];
d = pdist2(xtrue', a');

cvx_begin quiet
variable x(2);
minimize(0);
subject to
for i = 1:3
    norm(x - a(:, i)) <= d(i);
end
cvx_end
norm(x-xtrue)

%% SDP Relaxation
a = [1, -1, 0; 0, 0, 2];
lambda = rand(3,1);
lambda = lambda / sum(lambda);
xtrue = a * lambda;
d = pdist2(xtrue', a');

A= [1, 0, 1; 0, 1,1; 0,0,0] ;

cvx_begin sdp
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


%%
clear all
close all
a = [1, -1, 0; 0, 0, 2];
[X_1, Y_1] = meshgrid(linspace(-1.5, 1.5,50)',linspace(-.5,2.5,50)')
norm_SPD = zeros(size(X_1));
norm_reg= zeros(size(X_1));
%%
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
        
        norm_SPD(j,l) = norm(Z(end,1:2) - [X_1(j,l), Y_1(j,l)]);
        
        
        cvx_begin quiet
            variable x(2);
            minimize(0.0);
            subject to
            for i = 1:3
                norm(x - a(:, i),2) <= d(i);
            end
        cvx_end
        
        norm_reg(j, l) = norm(x-[X_1(j,l); Y_1(j,l)],2);
        
        

        
    end
    
end

%%
figure()
subplot(2,1,1)
mesh(X_1, Y_1, norm_SPD)
title('Magntitude of Error with SDP Relaxation')
xlabel('X Data')
ylabel('Y Data')
zlabel('Error')

subplot(2,1,2)
mesh(X_1, Y_1, norm_reg)
title('Magntitude of Error with SOCP Relaxation')
xlabel('X Data')
ylabel('Y Data')
zlabel('Error')




%%
% close all
% figure()
% subplot(2,1,1)
% plot3(X(:,1), X(:,2), norm_SPD)
% title('Magntitude of Error with SOCP Relaxation')
% hold on
% tmp = a
% tmp(:,4) = tmp(:,1)
% plot3(tmp(1,:), tmp(2,:), zeros(size(tmp(2,:))))
% hold off
% xlabel('X Data')
% ylabel('Y Data')
% zlabel('Error')
% 
% subplot(2,1,2)
% plot3(X(:,1), X(:,2), norm_reg)
% title('Magntitude of Error without SPD Relaxation')
% hold on
% 
% plot3(tmp(1,:), tmp(2,:), zeros(size(tmp(2,:))))
% xlabel('X Data')
% ylabel('Y Data')
% zlabel('Error')

%%
% a = [1, -1, 0; 0, 0, 2];
% norm_not_hull = zeros(size(nh_data,1),1);
% for j=1:size(nh_data, 1)
%     j
%     A= [1, 0, 1; 0, 1,1; 0,0,0] ;
%     d = pdist2(nh_data(j,:), a');
%     cvx_begin sdp quiet
%         variable Z(3,3) symmetric
%         minimize(0);
%         subject to
%
%             sum(dot(A(:,1)*A(:,1)', Z)) == 1;
%             sum(dot(A(:,2)*A(:,2)', Z)) == 1;
%             sum(dot(A(:,3)*A(:,3)', Z)) == 2;
%
%             for i = 1:3
%                 sum(dot([a(:, i); -1] * [a(:, i); -1]',  Z)) == d(i)^2;
%             end
%
%             Z >= 0;
%
%     cvx_end
%
%
%     norm_inhull(j) = norm(Z(end,1:2) - nh_data(j,:));
%
%
%
% end
%
% close all
% figure()
%
% plot3(nh_data(:,1), nh_data(:,2), norm_ih)
% title('Magntitude of Error for Points in the Hull')
% hold on
% tmp = a
% tmp(:,4) = tmp(:,1)
% plot3(tmp(1,:), tmp(2,:), zeros(size(a(2,:))))
% xlabel('X Data')
% ylabel('Y Data')
%
% %%
% close all
% figure()
% subplot(2,1,1)
% plot3(inhull(:,1), inhull(:,2), norm_inhull)

