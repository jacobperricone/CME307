function [ Z ] = CVXSolZ(a, A, d, idx)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
C = 10*rand(4,4);
cvx_expert true
cvx_begin sdp quiet
            variable Z(4,4) symmetric
            minimize(sum(dot(C, Z)) + idx(1)*trace(Z) + ...
                idx(2)* (sum(dot([a(:, 3); -1; 0]*[a(:, 3); -1; 0]', Z)) ...
                +sum(dot([a(:, 1); -1; 0]*[a(:, 1); -1; 0]', Z))))
       
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


end

