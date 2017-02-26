clear all
close all


%% Problem 5
a = [0, 1, 0; 0 0, 1];
b = [0,-1, 0; 0, 0,-1];
mu = [0, 10^-5];

x_initial = [.5;1.25]
x0_initial= .124;
ALPHA = .005;
MAX_ITER = 10000;
TOL = .0001;

f = @(x,x0, a, b) sum(log(1 + exp(-a'*x - x0))) + sum(log(1 + exp(b'*x + x0))) + mu(i) * norm(x)^2;
df = @(x,x0, a, b) ...
    [sum((-a(1,:)'.*exp(-a'*x - x0)) ./ (1 + exp(-a'*x - x0))) ...
    + sum((b(1,:)'.*exp(b'*x + x0)) ./ (1 + exp(b'*x + x0))); ... 
    sum((-a(2,:)'.*exp(-a'*x - x0)) ./ (1 + exp(-a'*x - x0))) ...
    + sum((b(2,:)'.*exp(b'*x + x0)) ./ (1 + exp(b'*x + x0))); ... a
    sum((-exp(-a'*x - x0)) ./ (1 + exp(-a'*x - x0))) ...
    + sum((exp(b'*x + x0)) ./ (1 + exp(b'*x + x0)))] + 2*mu(i)*[x;0];

df =@(x,x0,a,b) sum(reshape(sum(horzcat(-a(1,:)',b(1,:)', -a(2,:)',b(2,:)',ones(size(a',1),2)) ... 
.* horzcat(repmat(horzcat(exp(-a'*x - x0) ./ (1 + exp(-a'*x - x0)),exp(b'*x + x0) ./ (1 + exp(b'*x + x0))),1,2,1), ...
-exp(-a'*x - x0) ./ (1 + exp(-a'*x - x0)), exp(b'*x + x0) ./ (1 + exp(b'*x + x0)))),2,3))'


for i  = 1:2
%% Call SDM and ASDM
    [x_SDM, x0_SDM] = SDM(f, df, x_initial, x0_initial, a, b,ALPHA, MAX_ITER, TOL,1);
 
    [x_ASDM, x0_ASMD] = ASDM(f, df, x_initial, x0_initial, a, b,1/ALPHA, MAX_ITER, TOL,0);

    [x_CGD, x0_CGD] = CGD(f, df, x_initial, x0_initial, a, b, MAX_ITER, TOL,0);

    [x_CGD, x0_CGD] = BB1(f, df, x_initial, x0_initial, a, b, MAX_ITER, TOL,1,0);
end
