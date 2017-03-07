clear all
close all
%%
a = [1, -1, 0; 0, 0, 2];
idx =[1, 0; -1, 0; 0, 1; 0, -1];

DMethods_error_C = [];
DMethods_error_V = [];
for i = 1:3
%   
    Methods_error_C = [];
    Methods_error_V = [];
    [X_1, Y_1] = meshgrid(linspace(-1.5, 1.5,20)',linspace(-.5,2.5,20)');
    
    Cx1_error_SDP = zeros(size(X_1,1), size(Y_1,1));
    Cx2_error_SDP = zeros(size(X_1,1), size(Y_1,1));
    Ctotal_error_SDP = zeros(size(X_1,1), size(Y_1,1));


    Vx1_error_SDP = zeros(size(X_1,1), size(Y_1,1));
    Vx2_error_SDP = zeros(size(X_1,1), size(Y_1,1));
    Vtotal_error_SDP = zeros(size(X_1,1), size(Y_1,1));

    A1 = [1; 0; 0; 0];   A2 = [0; 1; 0; 0];   A3 = [1; 1; 0; 0];
    A = [A1, A2, A3];

    for j=1:size(X_1, 1)
            
        for l=1:size(Y_1,1)
        
            x1_center = [0,1];
            dist_x1_center = pdist2(x1_center, a(:,1:2)');
        
            x1_vertex = [.95,.05];
            dist_x1_vertex = pdist2(x1_vertex, a(:,1:2)');
        
            x2 = [X_1(j,l),Y_1(j,l)];
            dist_x2 = pdist2(x2, a(:,2:3)');
        
            dhat_center = pdist2(x1_center, x2);
            dhat_vertex = pdist2(x1_vertex, x2);
        
       
            d = [dist_x1_center(1), dist_x1_center(2), dist_x2(1), ...
            dist_x2(2), dhat_center];
        

            
            % Minimize the trace of Z
            Z = CVXSolZ(a, A, d, idx(i, :));

            x_1 = [Z(3,1), Z(3,2)]';
            x_2 = [Z(4,1), Z(4,2)]';




            Cx1_error_SDP(j,l) = norm(x_1 - x1_center);
            Cx2_error_SDP(j,l) = norm(x_2 - x2);
            Ctotal_error_SDP(j,l) = Cx1_error_SDP(j,l) +  Cx2_error_SDP(j,l);
            Cdistance(j,l) = dhat_center;
                
            Methods_error_C = [Methods_error_C, Cx1_error_SDP(j,l) +  Cx2_error_SDP(j,l)];    
                
            d = [dist_x1_vertex(1), dist_x1_vertex(2), dist_x2(1), ...
                dist_x2(2), dhat_vertex];
    
                
            % Minimize the trace of Z
            Z = CVXSolZ(a, A, d, idx(i, :));

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
        
            Vx1_error_SDP(j,l) = norm(x_1 - x1_vertex);
            Vx2_error_SDP(j,l) = norm(x_2 - x2);
            Vtotal_error_SDP(j,l) = Vx1_error_SDP(j,l) +  Vx2_error_SDP(j,l);
            Vdistance(j,l) = dhat_vertex;
            
            Methods_error_V = [Methods_error_V, Vx1_error_SDP(j,l) +  Vx2_error_SDP(j,l)];
        
        end
    end
 
    
    %%
    idx =[1, 0; -1, 0; 0, 1; 0, -1];
    col = idx(i,:);
    if col(1) == 1 || col(2) == 0
        addinfo = ' Objective is to Minimize C\cdot Z + tr(Z).';
    
    elseif col(1) == -1 col(2) == 0
        addinfo = ' Objective is to Minimize C\cdot Z - tr(Z).';
    
    elseif col(1) == 0 || col(2) == 1
        addinfo = ' Objective is to Minimize C\cdotZ + (d_{13} + d_{21})^2';
    
    else
        addinfo = ' Objective is to Minimize C\cdot Z - (d_{13} + d_{21})^2';
    end
    
    %%
    figure()
    subplot(3,1,1)

    mesh(X_1, Y_1, Cx1_error_SDP)
    colormap hsv
    alpha(.4)
    colorbar
    view(-30,30); camlight; axis image


    title(strcat('Magntitude of X_1 Error with X_1 Fixed at (x = 0, y = 1) in Hull.', addinfo), 'FontSize', 10)
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

    title(strcat('Magntitude of X_2 Error with X_1 Fixed at (x = 0, y = 1) in Hull', addinfo), 'FontSize', 10)
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



    title(strcat('Magntitude of Total Error with X_1 Fixed at (x = 0, y = 1) in Hull', addinfo), 'FontSize', 10)
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


    title(strcat('Magntitude of X_1 Error with X_1 Fixed at (x = .95, y = .05) in Hull', addinfo), 'FontSize', 10)
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


    title(strcat('Magntitude of X_2 Error with X_1 Fixed at (x = .95,y = .05 ) in Hull', addinfo), 'FontSize', 10)
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


    title(strcat('Magntitude of Total Error with X_1 Fixed at (x = .95,y = .05) in Hull', addinfo), 'FontSize', 10)
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
    
    
    DMethods_error_C = [DMethods_error_C; Methods_error_C];
    DMethods_error_V = [DMethods_error_V; Methods_error_V];

end
figure()
plot(1:length(DMethods_error_C(1, :)), DMethods_error_C(1, :))
hold on
plot(1:length(DMethods_error_C(2, :)), DMethods_error_C(2, :))
hold on
plot(1:length(DMethods_error_C(3, :)), DMethods_error_C(3, :))
title('Differences Across Methods in Estimating Fixed X = (0, 1)')
legend('With tr(Z)', 'With -tr(Z)', 'With (d_{13} + d_{21})^2')
hold off

figure()
plot(1:length(DMethods_error_V(1, :)), DMethods_error_V(1, :))
hold on
plot(1:length(DMethods_error_V(2, :)), DMethods_error_V(2, :))
hold on
plot(1:length(DMethods_error_V(3, :)), DMethods_error_V(3, :))
title('Differences Across Methods in Estimating Fixed X = (.95, .05)')
legend('With tr(Z)', 'With -tr(Z)', 'With (d_{13} + d_{21})^2')
