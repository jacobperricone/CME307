function [x] = SDMwFtwoXs(alpha, a, d, x1_0, x2_0)

x1(:, 1) = x1_0;
x2(:, 1) = x2_0;
%fun(1) = f0;
%out(1, :) =[x1(:, 1)', x2(:, 1)'];

for i = 1: 15000 -1 
    x1(:, i+1) = x1(:, i) - alpha * gradx1(a, d, x1(:, i), x2(:, i));
    x2(:, i+1) = x2(:, i) - alpha * gradx2(a, d, x1(:, i+1), x2(:, i));
    %fun(i+1)= fx1x2(x1(:, i+1), x2(:, i+1), a, d);
    %out(i+1, :) = [x1(:, i+1)', x2(:, i+1)'];
end

%x = [out(end, 1:4)];
x = [x1(:, end), x2(:, end)];
end
