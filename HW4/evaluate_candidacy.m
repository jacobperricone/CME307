function [out,points] = evaluate_candidacy(x, s, p, cutoff)
 
    out = p*diag(log(x'*s))' - sum(log(x.*s),1);
    points = x(:,out < cutoff);
    
    
end