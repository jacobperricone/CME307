clear all
close all


%% Problem 5
a = [0, 1, 0; 0 0, 1];
b = [0,-1, 0; 0, 0,-1];
mu = [0, 10^-5];

x_initial = [.5;1.25]
x0_initial= .124;
ALPHA = .005;
MAX_ITER = 25000;
TOL = .0001;

f = @(x,x0,  a, b, mu) sum(log(1 + exp(-a'*x - x0))) + sum(log(1 + exp(b'*x + x0))) + mu * norm(x)^2;
df = @(x,x0, a, b, mu) ...
    [sum((-a(1,:)'.*exp(-a'*x - x0)) ./ (1 + exp(-a'*x - x0))) ...
    + sum((b(1,:)'.*exp(b'*x + x0)) ./ (1 + exp(b'*x + x0))); ... 
    sum((-a(2,:)'.*exp(-a'*x - x0)) ./ (1 + exp(-a'*x - x0))) ...
    + sum((b(2,:)'.*exp(b'*x + x0)) ./ (1 + exp(b'*x + x0))); ... 
    sum((-exp(-a'*x - x0)) ./ (1 + exp(-a'*x - x0))) ...
    + sum((exp(b'*x + x0)) ./ (1 + exp(b'*x + x0)))] + 2*mu*[x;0];

%%
newa = vertcat(a,ones(1,size(a,2)))
tmpa = zeros(size(a,2), size(a,1)*(size(a,2) + 1));
tmpa(1,1:3) = newa(:,1)';
tmpa(2,4:6) = newa(:,2)'
tmpa(3,7:9) = newa(:,3)'
tmpa = reshape(newa*tmpa,3,3,size(a,2));

newb = vertcat(b,ones(1,size(b,2)))
tmpb = zeros(size(b,2), size(b,1)*(size(b,2) + 1));
tmpb(1,1:3) = newb(:,1)';
tmpb(2,4:6) = newb(:,2)'
tmpb(3,7:9) = newb(:,3)'
tmpb = reshape(newb*tmpb,3,3,size(b,2));


z = @(x,x0,a) (repmat(reshape(exp(-a'*x - x0)./( 1 + exp(-a'*x - x0)).^2,1,1,size(newa,2)),size(newa,1),size(newa,1),1));
zbar  = @(x,x0,b)(repmat(reshape(exp(b'*x + x0)./(1 + exp(b'*x + x0)).^2,1,1,size(newb,2)),size(newb,1),size(newb,1),1));


hessian = @(x,x0,a,b,mu)(sum(z(x,x0,a).*tmpa,3) + sum(zbar(x,x0,b).*tmpb,3))
%%

for i  = 1:2
% Call SDM and ASDM

    [x_SDM, x0_SDM] = SDM(f, df, x_initial, x0_initial, a, b,ALPHA, MAX_ITER, TOL,0,mu(i));
 
    [x_ASDM, x0_ASMD] = ASDM(f, df, x_initial, x0_initial, a, b,1/ALPHA, MAX_ITER, TOL,0,mu(i));

    [x_CGD, x0_CGD] = CGD(f, df, x_initial, x0_initial, a, b, MAX_ITER, TOL,0,mu(i));

    [x_CGD, x0_CGD] = BB1(f, df, x_initial, x0_initial, a, b, MAX_ITER, TOL,1,1,mu(i));
end