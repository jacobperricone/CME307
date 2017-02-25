clear all
close all


a = [1, -1, 0; 0, 0, 2];


x1 = [2*rand(25,1) + -1, 2*rand(25,1)];
x2 = [2*rand(25,1) + -1, 2*rand(25,1)];


dist_x1 = pdist2(x1, a(:,1:2)');
dist_x2 = pdist2(x2, a(:,2:3)');
dhat = pdist2(x1, x2);



%% Fix point in convex hull



%% Middle of HULL

clear all
close all


a = [1, -1, 0; 0, 0, 2];
%%

[X_1, Y_1] = meshgrid(linspace(-1.5, 1.5,35)',linspace(-.5,2.5,35)');

Cx1_error_SOCP = zeros(size(X_1,1), size(Y_1,1));
Cx2_error_SOCP = zeros(size(X_1,1), size(Y_1,1));
Ctotal_error_SOCP = zeros(size(X_1,1), size(Y_1,1));


Vx1_error_SOCP = zeros(size(X_1,1), size(Y_1,1));
Vx2_error_SOCP = zeros(size(X_1,1), size(Y_1,1));
Vtotal_error_SOCP = zeros(size(X_1,1), size(Y_1,1));



for j=1:size(X_1, 1)
    j;
    for l=1:size(Y_1,1)
        l;
        x1_center = [0,1];
        dist_x1_center = pdist2(x1_center, a(:,1:2)');
        
        x1_vertex = [.95,.05];
        dist_x1_vertex = pdist2(x1_vertex, a(:,1:2)');
        
        x2 = [X_1(j,l),Y_1(j,l)];
        dist_x2 = pdist2(x2, a(:,2:3)');
        
        dhat_center = pdist2(x1_center, x2);
        dhat_vertex = pdist2(x1_vertex, x2);
        
        cvx_begin quiet
            variables x_1(2) x_2(2)
            minimize 0
            subject to
            for i = 1:2
                norm(x_1 - a(:, i),2) <= dist_x1_center(i);
                norm(x_2 - a(:,i+1),2) <= dist_x2(i);
                norm(x_1 - x_2,2) <= dhat_center;
            end
        cvx_end
        
        
        
        Cx1_error_SOCP(j,l) = norm(x_1 - x1_center');
        Cx2_error_SOCP(j,l) = norm(x_2 - x2');
        Ctotal_error_SOCP(j,l) = Cx1_error_SOCP(j,l) +  Cx2_error_SOCP(j,l);
        Cdistance(j,l) = dhat_center;
        
        cvx_begin quiet
            variables x_1(2) x_2(2)
            minimize 0
            subject to
            for i = 1:2
                norm(x_1 - a(:, i),2) <= dist_x1_vertex(i);
                norm(x_2 - a(:,i+1),2) <= dist_x2(i);
                norm(x_1 - x_2,2) <= dhat_vertex;
            end
        cvx_end
        
        Vx1_error_SOCP(j,l) = norm(x_1 - x1_vertex');
        Vx2_error_SOCP(j,l) = norm(x_2 - x2');
        Vtotal_error_SOCP(j,l) = Vx1_error_SOCP(j,l) +  Vx2_error_SOCP(j,l);
        Vdistance(j,l) = dhat_vertex;
        
        
    end
end
 
%%
figure()
subplot(3,1,1)
surf(X_1, Y_1, Cx1_error_SOCP)
colormap hsv
alpha(.4)
colorbar
view(-30,30); camlight; axis image



title('Magntitude of X_1 Error with SOCP Relaxation X_1 Fixed at (x = 0, y = 1) in Hull', 'FontSize', 10)
xlabel('X Data')
ylabel('Y Data')
zlabel('Error')
[mx,k] = min(Cx1_error_SOCP(:));
[ix,jx] = ind2sub(size(Cx1_error_SOCP),k);
dim = [.10 .595 .3 .3];
str = strcat('Minimum Error of (', num2str(Cx1_error_SOCP(ix,jx)),')', ' at X_2: (', num2str(X_1(ix,jx)), ',', num2str(Y_1(ix,jx)),')');
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'FontSize',10);

hold on
tmp = a;
tmp(:,4) = tmp(:,1);
plot3(tmp(1,:), tmp(2,:), zeros(size(tmp(2,:))), '-r')
plot3(x1_center(1), x1_center(2), 0, '-*b','MarkerSize',10)
xlabel('X Data')
ylabel('Y Data')
hold off


subplot(3,1,2)
surf(X_1, Y_1, Cx2_error_SOCP)
colormap hsv
alpha(.4)
colorbar
view(-30,30); camlight; axis image


title('Magntitude of X_2 Error with SOCP Relaxation X_1 Fixed at (x = 0, y = 1) in Hull', 'FontSize', 10)
xlabel('X Data')
ylabel('Y Data')
zlabel('Error')

[mx,k] = min(Cx2_error_SOCP(:));
[ix,jx] = ind2sub(size(Cx2_error_SOCP),k);
dim = [.10 .295 .3 .3];
str = strcat('Minimum Error of (', num2str(Cx2_error_SOCP(ix,jx)),')', ' at X_2: (', num2str(X_1(ix,jx)), ',', num2str(Y_1(ix,jx)),')');
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'FontSize',10);

hold on
tmp = a;
tmp(:,4) = tmp(:,1);
plot3(tmp(1,:), tmp(2,:), zeros(size(tmp(2,:))), '-r')
plot3(x1_center(1), x1_center(2), 0, '-*b','MarkerSize',10)
xlabel('X Data')
ylabel('Y Data')
hold off


subplot(3,1,3)
surf(X_1, Y_1, Ctotal_error_SOCP)
colormap hsv
alpha(.4)
colorbar
view(-30,30); camlight; axis image



title('Magntitude of Total Error with SOCP Relaxation X_1 Fixed at (x = 0, y = 1) in Hull', 'FontSize', 10)
xlabel('X Data')
ylabel('Y Data')
zlabel('Error')

[mx,k] = min(Ctotal_error_SOCP(:));
[ix,jx] = ind2sub(size(Ctotal_error_SOCP),k);
dim = [.10 .009 .3 .3];
str = strcat('Minimum Error of (', num2str(Ctotal_error_SOCP(ix,jx)),')', ' at X_2: (', num2str(X_1(ix,jx)), ',', num2str(Y_1(ix,jx)),')');
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'FontSize',10);

hold on
tmp = a;
tmp(:,4) = tmp(:,1);
plot3(tmp(1,:), tmp(2,:), zeros(size(tmp(2,:))), '-r')
plot3(x1_center(1), x1_center(2), 0, '-*b','MarkerSize',10)
xlabel('X Data')
ylabel('Y Data')
hold off


%
figure()
subplot(3,1,1)
surf(X_1, Y_1, Vx1_error_SOCP)
colormap hsv
alpha(.4)
colorbar
view(-30,30); camlight; axis image




title('Magntitude of X_1 Error with SOCP Relaxation X_1 Fixed at (x = .95, y = .05) in Hull', 'FontSize', 10)
xlabel('X Data')
ylabel('Y Data')
zlabel('Error')



[mx,k] = min(Vx1_error_SOCP(:));
[ix,jx] = ind2sub(size(Vx1_error_SOCP),k);
dim = [.10 .595 .3 .3];
str = strcat('Minimum Error of (', num2str(Vx1_error_SOCP(ix,jx)),')', ' at X_2: (', num2str(X_1(ix,jx)), ',', num2str(Y_1(ix,jx)),')');
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'FontSize',10);

hold on
tmp = a;
tmp(:,4) = tmp(:,1);
plot3(tmp(1,:), tmp(2,:), zeros(size(tmp(2,:))), '-r')
plot3(x1_vertex(1), x1_vertex(2), 0, '-*b','MarkerSize',10)
xlabel('X Data')
ylabel('Y Data')
hold off


subplot(3,1,2)


surf(X_1, Y_1, Vx2_error_SOCP)
colormap hsv
alpha(.4)
colorbar
view(-30,30); camlight; axis image


title('Magntitude of X_2 Error with SOCP Relaxation X_1 Fixed at (x = .95,y = .05 ) in Hull', 'FontSize', 10)
xlabel('X Data')
ylabel('Y Data')
zlabel('Error')
[mx,k] = min(Vx2_error_SOCP(:));
[ix,jx] = ind2sub(size(Vx2_error_SOCP),k);
dim = [.10 .295 .3 .3];
str = strcat('Minimum Error of (', num2str(Vx2_error_SOCP(ix,jx)),')', ' at X_2: (', num2str(X_1(ix,jx)), ',', num2str(Y_1(ix,jx)),')');
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'FontSize',10);

hold on
tmp = a;
tmp(:,4) = tmp(:,1);
plot3(tmp(1,:), tmp(2,:), zeros(size(tmp(2,:))), '-r')
plot3(x1_vertex(1), x1_vertex(2), 0, '-b','MarkerSize',10)
xlabel('X Data')
ylabel('Y Data')
hold off


subplot(3,1,3)
surf(X_1, Y_1, Vtotal_error_SOCP)
colormap hsv
alpha(.4)
colorbar
view(-30,30); camlight; axis image


title('Magntitude of Total Error with SOCP Relaxation X_1 Fixed at (x = .95,y = .05) in Hull', 'FontSize', 10)
xlabel('X Data')
ylabel('Y Data')
zlabel('Error')
[mx,k] = min(Vtotal_error_SOCP(:));
[ix,jx] = ind2sub(size(Vtotal_error_SOCP),k);
dim = [.1 .011 .3 .3];
str = strcat('Minimum Error of (', num2str(Vtotal_error_SOCP(ix,jx)),')', ' at X_2: (', num2str(X_1(ix,jx)), ',', num2str(Y_1(ix,jx)),')');
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'FontSize',10);


hold on
tmp = a;
tmp(:,4) = tmp(:,1);
plot3(tmp(1,:), tmp(2,:), zeros(size(tmp(2,:))), '-r')
plot3(x1_vertex(1), x1_vertex(2), 0, '-*b','MarkerSize',10)
xlabel('X Data')
ylabel('Y Data')
hold off


%%

clear all

a = [1, -1, 0; 0, 0, 2];



%
[X_1, Y_1] = meshgrid(linspace(-1.5, 1.5,35)',linspace(-.5,2.5,35)');

Cx1_error_SDP = zeros(size(X_1,1), size(Y_1,1));
Cx2_error_SDP = zeros(size(X_1,1), size(Y_1,1));
Ctotal_error_SDP = zeros(size(X_1,1), size(Y_1,1));


Vx1_error_SDP = zeros(size(X_1,1), size(Y_1,1));
Vx2_error_SDP = zeros(size(X_1,1), size(Y_1,1));
Vtotal_error_SDP = zeros(size(X_1,1), size(Y_1,1));

A1 = [1; 0; 0; 0];   A2 = [0; 1; 0; 0];   A3 = [1; 1; 0; 0];
A = [A1, A2, A3];
a = [1, -1, 0; 0, 0, 2];
for j=1:size(X_1, 1)
    j;
    for l=1:size(Y_1,1)
        
        x1_center = [0,1];
        dist_x1_center = pdist2(x1_center, a(:,1:2)');
        
        x1_vertex = [.95,.05];
        dist_x1_vertex = pdist2(x1_vertex, a(:,1:2)');
        
        x2 = [X_1(j,l),Y_1(j,l)];
        dist_x2 = pdist2(x2, a(:,2:3)');
        
        dhat_center = pdist2(x1_center, x2);
        dhat_vertex = pdist2(x1_vertex, x2);
        
        cvx_begin sdp quiet
            variable Z(4,4) symmetric
            minimize(0);
            subject to

            sum(dot(A(:,1)*A(:,1)', Z)) == 1;
            sum(dot(A(:,2)*A(:,2)', Z)) == 1;
            sum(dot(A(:,3)*A(:,3)', Z)) == 2;
            
            for i = 1:2
                
                sum(dot([a(:, i); -1; 0] * [a(:, i); -1; 0]',  Z))...
                    == dist_x1_center(i)^2;
                sum(dot([a(:, i+1); 0;  -1] * [a(:, i+1); 0; -1]', Z))...
                    == dist_x2(i)^2;
            end

            sum(dot([0; 0; 1; -1] * [0; 0; 1; -1]', Z)) == dhat_center^2 ;

            Z >= 0;

        cvx_end
        
        x_1 = [Z(3,1), Z(3,2)]';
        x_2 = [Z(4,1), Z(4,2)]';
        
        

        
        
            
        Cx1_error_SDP(j,l) = norm(x_1 - x1_center');
        Cx2_error_SDP(j,l) = norm(x_2 - x2');
        Ctotal_error_SDP(j,l) = Cx1_error_SDP(j,l) +  Cx2_error_SDP(j,l);
        Cdistance(j,l) = dhat_center;
        

        cvx_begin sdp quiet
            variable Z(4,4) symmetric
            minimize(0);
            subject to

            sum(dot(A(:,1)*A(:,1)', Z)) == 1;
            sum(dot(A(:,2)*A(:,2)', Z)) == 1;
            sum(dot(A(:,3)*A(:,3)', Z)) == 2;

            for i = 1:2
                sum(dot([a(:, i); -1; 0] * [a(:, i); -1; 0]',  Z))...
                    == dist_x1_vertex(i)^2;
                sum(dot([a(:, i+1); 0;  -1] * [a(:, i+1); 0; -1]', Z))...
                    == dist_x2(i)^2;
            end

            sum(dot([0; 0; 1; -1] * [0; 0; 1; -1]', Z)) == dhat_vertex^2 ;

            Z >= 0;

        cvx_end
        
        
        x_1 = [Z(3,1), Z(3,2)]';
        x_2 = [Z(4,1), Z(4,2)]';
        

        
        if isnan(x_1(1)) || isnan(x_1(2))
            disp('YOOx1')
            x_1;
        end
        
        if isnan(x_1(1)) || isnan(x_1(2))
            disp('YOOx2')
            x_2;
        end
        
        Vx1_error_SDP(j,l) = norm(x_1 - x1_vertex');
        Vx2_error_SDP(j,l) = norm(x_2 - x2');
        Vtotal_error_SDP(j,l) = Vx1_error_SDP(j,l) +  Vx2_error_SDP(j,l);
        Vdistance(j,l) = dhat_vertex;
        
    end
end
 
%%
figure()
subplot(3,1,1)

surf(X_1, Y_1, Cx1_error_SDP)
colormap hsv
alpha(.4)
colorbar
view(-30,30); camlight; axis image


title('Magntitude of X_1 Error with SDP Relaxation X_1 Fixed at (x = 0, y = 1) in Hull', 'FontSize', 10)
xlabel('X Data')
ylabel('Y Data')
zlabel('Error')

[mx,k] = min(Cx1_error_SDP(:));
[ix,jx] = ind2sub(size(Cx1_error_SDP),k);
dim = [.10 .605 .3 .3];
str = strcat('Minimum Error of (', num2str(Cx1_error_SDP(ix,jx)),')', ' at X_2: (', num2str(X_1(ix,jx)), ',', num2str(Y_1(ix,jx)),')');
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'FontSize',10);

      
 

hold on
tmp = a;
tmp(:,4) = tmp(:,1);
plot3(tmp(1,:), tmp(2,:), zeros(size(tmp(2,:))), '-r')
plot3(x1_center(1), x1_center(2), 0, '-*b','MarkerSize',10)
xlabel('X Data')
ylabel('Y Data')
hold off


subplot(3,1,2)
surf(X_1, Y_1, Cx2_error_SDP)
colormap hsv
alpha(.4)
colorbar
view(-30,30); camlight; axis image

title('Magntitude of X_2 Error with SDP Relaxation X_1 Fixed at (x = 0, y = 1) in Hull', 'FontSize', 10)
xlabel('X Data')
ylabel('Y Data')
zlabel('Error')

[mx,k] = min(Cx2_error_SDP(:));
[ix,jx] = ind2sub(size(Cx2_error_SDP),k);
dim = [.10 .295 .3 .3];
str = strcat('Minimum Error of (', num2str(Cx2_error_SDP(ix,jx)),')', ' at X_2: (', num2str(X_1(ix,jx)), ',', num2str(Y_1(ix,jx)),')');
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'FontSize',10);

hold on
tmp = a;
tmp(:,4) = tmp(:,1);
plot3(tmp(1,:), tmp(2,:), zeros(size(tmp(2,:))), '-r')
plot3(x1_center(1), x1_center(2), 0, '-*b','MarkerSize',10)
xlabel('X Data')
ylabel('Y Data')
hold off



subplot(3,1,3)
surf(X_1, Y_1, Ctotal_error_SDP)
colormap hsv
alpha(.4)
colorbar
view(-30,30); camlight; axis image



title('Magntitude of Total Error with SDP Relaxation X_1 Fixed at (x = 0, y = 1) in Hull', 'FontSize', 10)
xlabel('X Data')
ylabel('Y Data')
zlabel('Error')

[mx,k] = min(Ctotal_error_SDP(:));
[ix,jx] = ind2sub(size(Ctotal_error_SDP),k);
dim = [.10 .009 .3 .3];
str = strcat('Minimum Error of (', num2str(Ctotal_error_SDP(ix,jx)),')', ' at X_2: (', num2str(X_1(ix,jx)), ',', num2str(Y_1(ix,jx)),')');
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'FontSize',10);

hold on
tmp = a;
tmp(:,4) = tmp(:,1);
plot3(tmp(1,:), tmp(2,:), zeros(size(tmp(2,:))), '-r')
plot3(x1_center(1), x1_center(2), 0, '-*b','MarkerSize',10)
xlabel('X Data')
ylabel('Y Data')
hold off
%
figure()
subplot(3,1,1)
surf(X_1, Y_1, Vx1_error_SDP)
colormap hsv
alpha(.4)
colorbar
view(-30,30); camlight; axis image


title('Magntitude of X_1 Error with SDP Relaxation X_1 Fixed at (x = .95, y = .05) in Hull', 'FontSize', 10)
xlabel('X Data')
ylabel('Y Data')
zlabel('Error')


[mx,k] = min(Vx1_error_SDP(:));
[ix,jx] = ind2sub(size(Vx1_error_SDP),k);
dim = [.10 .595 .3 .3];
str = strcat('Minimum Error of (', num2str(Vx1_error_SDP(ix,jx)),')', ' at X_2: (', num2str(X_1(ix,jx)), ',', num2str(Y_1(ix,jx)),')');
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'FontSize',10);

hold on
tmp = a;
tmp(:,4) = tmp(:,1);
plot3(tmp(1,:), tmp(2,:), zeros(size(tmp(2,:))), '-r')
plot3(x1_vertex(1), x1_vertex(2), 0, '-*b','MarkerSize',10)
xlabel('X Data')
ylabel('Y Data')
hold off



subplot(3,1,2)
surf(X_1, Y_1, Vx2_error_SDP)
colormap hsv
alpha(.4)
colorbar
view(-30,30); camlight; axis image


title('Magntitude of X_2 Error with SDP Relaxation X_1 Fixed at (x = .95,y = .05 ) in Hull', 'FontSize', 10)
xlabel('X Data')
ylabel('Y Data')
zlabel('Error')


[mx,k] = min(Vx2_error_SDP(:));
[ix,jx] = ind2sub(size(Vx2_error_SDP),k);
dim = [.10 .295 .3 .3];
str = strcat('Minimum Error of (', num2str(Vx2_error_SDP(ix,jx)),')', ' at X_2: (', num2str(X_1(ix,jx)), ',', num2str(Y_1(ix,jx)),')');
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'FontSize',10);

hold on
tmp = a;
tmp(:,4) = tmp(:,1);
plot3(tmp(1,:), tmp(2,:), zeros(size(tmp(2,:))), '-r')
plot3(x1_vertex(1), x1_vertex(2), 0, '-*b','MarkerSize',10)
xlabel('X Data')
ylabel('Y Data')
hold off


subplot(3,1,3)
surf(X_1, Y_1, Vtotal_error_SDP)
colormap hsv
alpha(.4)
colorbar
view(-30,30); camlight; axis image


title('Magntitude of Total Error with SDP Relaxation X_1 Fixed at (x = .95,y = .05) in Hull', 'FontSize', 10)
xlabel('X Data')
ylabel('Y Data')
zlabel('Error')
[mx,k] = min(Vtotal_error_SDP(:));
[ix,jx] = ind2sub(size(Vtotal_error_SDP),k);
dim = [.10 .011 .3 .3];
str = strcat('Minimum Error of (', num2str(Vtotal_error_SDP(ix,jx)),')', ' at X_2: (', num2str(X_1(ix,jx)), ',', num2str(Y_1(ix,jx)),')');
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'FontSize',10);


hold on
tmp = a;
tmp(:,4) = tmp(:,1);
plot3(tmp(1,:), tmp(2,:), zeros(size(tmp(2,:))), '-r')
plot3(x1_vertex(1), x1_vertex(2), 0, '-*b','MarkerSize',10)
xlabel('X Data')
ylabel('Y Data')
hold off