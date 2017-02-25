clear all
close all



w = sym('w',[1,3 ])
what = sym('what',[1,3])
x = sym('x', [1,2])
d = sym('x', [1,3])
A = sym('a_%d_%d', [3,3])


S = solve(w(1) + w(3) + sum((A(1,:).^2).*what) == -x(1)^2 ... 
     , w(3) + sum((A(1,:).*A(2,:)).*what) == -x(2)*x(1) ... 
     , sum(A(1,:).*what) == -x(1) ... 
     , sum(A(2,:).*what) == -x(1) ...
     , w(2) + w(3) + sum((A(2,:).^2).*what) == -x(2)^2 ...
     , sum(what) == -1 , [w, what])
 
 
 
 
 
 
 syms x1 y1 x2 y2 
 
 
 S = solve(x1 + y1 == 1, 3*x2 + y2 == 1, x1 + x2 == 1, y1 + y2 == 1, [x1, y1, x2, y2])
 
 
 X = linsolve([1, 1, 0, 0; 0, 0, 3, 1; 1,0, 1, 0; 0, 1, 0, 1], [1;1;1;1])