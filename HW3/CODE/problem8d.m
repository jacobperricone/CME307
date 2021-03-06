%%
% Author: Jacob Perricone
%
%
clear all
close all
%% HW 1 AFFINE TRANFORM


% set a points
a = [1, -1, 0; 0, 0, 2];

% create matrix
A = zeros(3,3,3);
% fill three-d matrix
for i=1:3
    A(:,:,i) = [a(:,i); -1]*[a(:,i); -1]';
end


%% make x
X = eye(3,3)
% set y
X(end,end) = 2;
% set initial points
x_initial = [-.5;1]
X(1:2,end) = x_initial
X(end,1:2)= x_initial'

% save true
x_true = [.2;.5]
% calculate distance
b = (pdist2(x_true', a').^2)';
% have mus
mus = linspace(0,10e-4,5)
% define functiona nd gradient for two versions
f = @(A,X,mu)(.5*norm(squeeze(sum(dot(A, repmat(X,1,1,size(A,3))))) - b,2)^2 - mu*log(det(X)))
df_affine = @(A,X,mu)(X*sum(A .* repmat(reshape(squeeze(sum(dot(A,repmat(X,1,1, size(A,3))))) - b,1,1,size(A,3)),3,3,1),3)*X - mu*X)
df_reg = @(A,X,mu)(sum(A .* repmat(reshape(squeeze(sum(dot(A,repmat(X,1,1, size(A,3))))) - b,1,1,size(A,3)),3,3,1),3) - mu*inv(X))




%%

% Define alpha
ALPHA = .01;
MAX_ITER = 50000;
TOL = 1.e-10;


%%
X_final_reg = zeros(size(mus,2),2,1)
X_final_affine = zeros(size(mus,2),2,1)
for i=1:size(mus,2)
    [tmp, tmp1] = SDP_REGRESSION(f, df_reg, A,X,mus(i),x_true,TOL, 500000, 1/ALPHA,a,0, 'No Scaling');
    X_final_reg(i,:,1) = tmp
    [tmp, tmp1] = SDP_REGRESSION(f, df_affine, A,X,mus(i),x_true,TOL, MAX_ITER, 1/ALPHA,a,0, 'Affine Scaling');
    X_final_affine(i,:,1) = tmp
end

%% HW 2 AFFINE TRANSFOR
% A1 = [1; 0; 0; 0];   A2 = [0; 1; 0; 0];   A3 = [1; 1; 0; 0];
% A = [A1, A2, A3];
a = [1, -1, 0; 0, 0, 2];
A = zeros(4,4,4)
for i=1:5
    if i<=2
        A(:,:,i) = [a(:,i);-1;0]*[a(:,i);-1;0]';
    elseif i== 3 || i == 4
        A(:,:,i) = [a(:,i-1);0;-1]*[a(:,i-1);0;-1]';
    else
        A(:,:,i) = [0;0;1;-1]*[0;0;1;-1]';
    end
    
    
    
    
end
%%
    
x1_true = [.8;.1];
x2_true = [0; 1.8];

x1_initial = [.8; .3];
x2_initial = [0; 1.5];

b(1:2) = (pdist2(x1_true', a(:,1:2)')').^2;
b(3:4) = (pdist2(x2_true', a(:,2:3)')').^2;
b(5) = (pdist2(x2_true', x1_true')').^2;

X = eye(4,4);
X(end-1:end,1:2) = [x1_initial,x2_initial]';
X(1:2, end-1:end) = [x1_initial,x2_initial];
X(end,end) = 3.5

f = @(A,X,mu)(.5*norm(squeeze(sum(dot(A, repmat(X,1,1,size(A,3))))) - b,2)^2 - mu*log(det(X)))
df_affine = @(A,X,mu)(X*sum(A .* repmat(reshape(squeeze(sum(dot(A,repmat(X,1,1, size(A,3))))) - b,1,1,size(A,3)),4,4,1),3)*X - mu*X)
df_reg = @(A,X,mu)(sum(A .* repmat(reshape(squeeze(sum(dot(A,repmat(X,1,1, size(A,3))))) - b,1,1,size(A,3)),4,4,1),3) - mu*inv(X))

mus = linspace(0,10e-4,2)
%%
ALPHA = .05
MAX_ITER = 100000;
TOL = 1.e-6;
x_true = horzcat(x1_true,x2_true);
%%
X_final = zeros(size(mus,2),2,1)
for i=1:size(mus,2)
[tmp1, tmp2,tmp3] = SDP_REGRESSION_HW2(f, df_reg, A,X,mus(i),x_true,TOL, MAX_ITER, 1/ALPHA,a,0, 'No Scaling')
[tmp1, tmp2,tmp3] = SDP_REGRESSION_HW2(f, df_affine, A,X,mus(i),x_true,TOL, MAX_ITER, 1/ALPHA,a,1, 'Affine Scaling');

end

%%



