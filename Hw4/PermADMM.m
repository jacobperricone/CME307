function [ ] = PermADMM(idx, beta, iter, A)

x1(1) = rand(1);
x2(1) = rand(1);
x3(1) = rand(1);
y(:, 1) = rand(3,1);
syms x_1 x_2 x_3

for k=1:iter
    order = randperm(3,3);
        
    if order(1)==1
        x1(k+1) = solve(idx*x_1 + A(:, 1)'*y(:, k) + beta * A(:, 1)'*(A(:, 1)*x_1 + ...
                A(:, 2)*x2(k) + A(:, 3)*x3(k)) == 0, x_1);
        if order(2)==2
            x2(k+1) = solve(idx*x_2+A(:, 2)'*y(:, k) + beta * A(:, 2)'*(A(:, 1)*x1(k+1) + ...
                    A(:, 2)*x_2 + A(:, 3)*x3(k)) == 0, x_2);

            x3(k+1) = solve(idx*x_3+A(:, 3)'*y(:, k) + beta * A(:, 3)'*(A(:, 1)*x1(k+1) + ...
                A(:, 2)*x2(k+1) + A(:, 3)*x_3) == 0, x_3);

        else
            x3(k+1) = solve(idx*x_3+A(:, 3)'*y(:, k) + beta * A(:, 3)'*(A(:, 1)*x1(k+1) + ...
                A(:, 2)*x2(k) + A(:, 3)*x_3) == 0, x_3);
                
            x2(k+1) = solve(idx*x_2+A(:, 2)'*y(:, k) + beta * A(:, 2)'*(A(:, 1)*x1(k+1) + ...
                    A(:, 2)*x_2 + A(:, 3)*x3(k+1)) == 0, x_2);
  
        end
        
        
    elseif order(1) == 2
            
        x2(k+1) = solve(idx*x_2+A(:, 2)'*y(:, k) + beta * A(:, 2)'*(A(:, 1)*x1(k) + ...
                    A(:, 2)*x_2 + A(:, 3)*x3(k)) == 0, x_2);
                
               
        if order(2) == 1
            x1(k+1) = solve(idx*x_1+A(:, 1)'*y(:, k) + beta * A(:, 1)'*(A(:, 1)*x_1 + ...
                    A(:, 2)*x2(k+1) + A(:, 3)*x3(k)) == 0, x_1);
                

            x3(k+1) = solve(idx*x_3+A(:, 3)'*y(:, k) + beta * A(:, 3)'*(A(:, 1)*x1(k+1) + ...
                    A(:, 2)*x2(k+1) + A(:, 3)*x_3) == 0, x_3);

               

        else
            x3(k+1) = solve(idx*x_3+A(:, 3)'*y(:, k) + beta * A(:, 3)'*(A(:, 1)*x1(k) + ...
                    A(:, 2)*x2(k+1) + A(:, 3)*x_3) == 0, x_3);
                
            x1(k+1) = solve(idx*x_1+A(:, 1)'*y(:, k) + beta * A(:, 1)'*(A(:, 1)*x_1 + ...
                    A(:, 2)*x2(k+1) + A(:, 3)*x3(k+1)) == 0, x_1);
  
        end
        
    else
            
        x3(k+1) = solve(idx*x_3+A(:, 3)'*y(:, k) + beta * A(:, 3)'*(A(:, 1)*x1(k) + ...
                    A(:, 2)*x2(k) + A(:, 3)*x_3) == 0, x_3);
            
        if order(2) == 1
            x1(k+1) = solve(idx*x_1+A(:, 1)'*y(:, k) + beta * A(:, 1)'*(A(:, 1)*x_1 + ...
                        A(:, 2)*x2(k) + A(:, 3)*x3(k+1)) == 0, x_1);
                    
                    x2(k+1) = solve(idx*x_2+A(:, 2)'*y(:, k) + beta * A(:, 2)'*(A(:, 1)*x1(k+1) + ...
                        A(:, 2)*x_2 + A(:, 3)*x3(k+1)) == 0, x_2);
                
        else
                    
            x2(k+1) = solve(idx*x_2+A(:, 2)'*y(:, k) + beta * A(:, 2)'*(A(:, 1)*x1(k) + ...
                        A(:, 2)*x_2 + A(:, 3)*x3(k+1)) == 0, x_2);
                    
            x1(k+1) = solve(idx*x_1+A(:, 1)'*y(:, k) + beta * A(:, 1)'*(A(:, 1)*x_1 + ...
                        A(:, 2)*x2(k+1) + A(:, 3)*x3(k+1)) == 0, x_1);
        end
    end
  
    y(:, k+1)=y(:, k)+beta * (A(:, 1)*x1(k+1) + A(:, 2)*x2(k+1) + A(:, 3)*x3(k+1));
            
            
            
            
      
    if mod(k,10) == 0
        disp(sprintf('----- beta: %d  ------', beta));
        disp(sprintf('------- Iteration: %d ----------', k));
        disp('      x1        x2        x3')
        disp([  x1(k), x2(k), x3(k)])
            
    end

end

if idx == 0
    add = ' without objective function';
else
    add = ' with objective function';
end

figure()
subplot(3,1,1)
plot(1:k+1, x1, '*', 1:k+1, zeros(1, k+1), '-')
title(strcat('Permuation based estimate of x_1 with beta=', num2str(beta), add))

subplot(3,1,2)
plot(1:k+1, x2, '*', 1:k+1, zeros(1, k+1), '-')
title(strcat('Permuation based estimate of x_1 with beta=', num2str(beta), add))

subplot(3,1,3)
plot(1:k+1, x3, '*', 1:k+1, zeros(1, k+1), '-')
title(strcat('Permuation based estimate of x_1 with beta=', num2str(beta), add))
hold off