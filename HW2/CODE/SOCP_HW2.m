clear all
close all

% Problem 6(b)
mu = [0, 10^-5];
a = [0, 1, 0; 0 0, 1];
b = [0,-1, 0; 0, 0, -1];


%%
% Assumption for x0
for j = 1:2
    cvx_begin quiet
        variable B
        variable x1 
        variable x2
        variable x0
        opt = B + mu(j) .* dot(x,x);

        minimize(opt);
        subject to

        B >= 0;


        x0 + B >=1
        x2 + x0 + B >= 1
        x2 + x0 + B >=1 
        
        x0 - B <= -1
        -x1 + x0 - B <= -1
        -x2 + x0 - B <= -1

        B >= 0
    cvx_end
    
    ex1{j} = x1;
    ex2{j} = x1;
    ex0{j} = x0;
    
    eB{j} = B;
end




% 
% tmp = {};
% for n=1:5
%     
%     for j = 1:2
%         cvx_begin sdp quiet
%         variable B
%         variable x(2)
%         variable x0
%         opt = B + mu(j) .* dot(x,x);
%         
%         minimize(opt);
%         subject to
%         
%         B >= 0;
%         
%         for i = 1:3
%             dot(a(:, i), x) + x0 + B >= 1;
%             dot(b(:, i), x) + x0 - B <= -1;
%         end
%         
%         cvx_end
%         
%         S.x{j} = x;
%         
%         S.opt{j} = opt;
%         S.B{j} = B
%     end
%     tmp{n}  =S;
%     
% end
%%


