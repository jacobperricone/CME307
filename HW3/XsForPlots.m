function [x] = XsForPlots( idx, a, d, A, ALPHA)

if idx == 1
    x_1 = [10*rand(1) - 1, 10*rand(1) - 1]';
    x_2 = [10*rand(1) - 1, 10*rand(1) - 1]';

elseif idx == 2
    
   cvx_begin sdp quiet
            variable Z(4,4) symmetric
            minimize(0);
            subject to

            sum(dot(A(:,1)*A(:,1)', Z)) == 1;
            sum(dot(A(:,2)*A(:,2)', Z)) == 1;
            sum(dot(A(:,3)*A(:,3)', Z)) == 2;
            
            for i = 1:2
                
                sum(dot([a(:, i); -1; 0] * [a(:, i); -1; 0]',  Z))...
                    == d(i)^2;
                sum(dot([a(:, i+1); 0;  -1] * [a(:, i+1); 0; -1]', Z))...
                    == d(i+2)^2;
            end

            sum(dot([0; 0; 1; -1] * [0; 0; 1; -1]', Z)) == d(5)^2 ;

            Z >= 0;

        cvx_end
        
        x_1 = [Z(3,1), Z(3,2)]';
        x_2 = [Z(4,1), Z(4,2)]';
  
end

x = SDMwFtwoXs(ALPHA, a, d, x_1, x_2);


end

