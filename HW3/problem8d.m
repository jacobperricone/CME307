%%
clear all
close all
%% HW 1


a = [1, -1, 0; 0, 0, 2];

A = zeros(3,3,3);

for i=1:3
    A(:,:,i) = [a(:,i); -1]*[a(:,i); -1]'
end

X = eye(3,3)
X(end,end) = 2;


x_initial = [-.5;1]
X(1:2,end) = x_initial
X(end,1:2)= x_initial'

x_true = [.2;.5]
b = (pdist2(x_true', a').^2)';
mus = linspace(0,10e-4,5)



f = @(A,X,mu)(.5*norm(squeeze(sum(dot(A, repmat(X,1,1,size(A,3))))) - b,2)^2 - mu*log(det(X)))
df = @(A,X,mu)(X*sum(A .* repmat(reshape(squeeze(sum(dot(A,repmat(X,1,1, 3)))) - b,1,1,3),3,3,1),3)*X - mu*X)




%%


ALPHA = .001;
MAX_ITER = 15000;
TOL = 1.e-10;


%%
X_final = zeros(size(mus,2),2,1)
for i=1:size(mus,2)
    tmp = SDP_REGRESSION(f, df, A,X,mus(i),x_true,TOL, MAX_ITER, 1/ALPHA,a,1);
    X_final(i,:,1) = tmp
end




