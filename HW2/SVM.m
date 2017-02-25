clear all
close all

% Problem 6(b)
mu = [0, 10^-5, 1, 2, 10, 100,100];
a = [0, 1, 0; 0 0, 1];
b = [0,-1, 0; 0, 0,-1];


%% Find optimal hyperplane for different Beta
for j = 1:size(mu,2)
    cvx_begin 
    variable B
    variable x1
    variable x2
    variable x0
    mu(j)
    opt = B + mu(j)*power(x1,2) + mu(j)*power(x2,2);
    
    minimize(opt);
    subject to
    
    B >= 0;
    
    
    x0 + B >=1
    x1 + x0 + B >= 1
    x2 + x0 + B >=1
    
    x0 - B <= -1
    -x1 + x0 - B <= -1
    -x2 + x0 - B <= -1
    
    B >= 0
    cvx_end
    
    ex1{j} = x1;
    ex2{j} = x2;
    ex{j} = [x1, x2];
    ex0{j} = x0;
    
    eB{j} = B;
end
%% Eliminate the Beta to See if results Change

cvx_begin quiet
    variable B
    variable x1
    variable x2
    variable x0
    opt = B
    minimize(opt);
    subject to

    B >= 0;
    x0 + B >=1
    x2 + x0 + B >= 1
    x1 + x0 + B >=1

    x0 - B <= -1
    -x1 + x0 - B <= -1
    -x2 + x0 - B <= -1
cvx_end



%%

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


