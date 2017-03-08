function alpha_new = choose_alpha(alpha, d,x, x0, a,b,f)
% finds appropriate alpha for next step
while f(x, x0, alpha, b) <= f(x + alpha * d(1:end-1), x0 + alpha*d(end), a, b) 
    alpha = alpha / 2;
 
end
alpha_new = alpha; 

end