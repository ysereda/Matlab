function x = Newton(fhandle,x,maxits,gtol,xtol,xstar)
% Computes the minimizer of f(x) using Newton's method
% Input:
%    fhandle - function handle for the function, gradient and Hessian
%    evaluation
%    x - initial guess
%    maxits - maximum number of iterations
%    gtol - tolerance for gradient
%    xtol - tolerance for step
%    xstar - true minimizer
% Output:
%    x - computed minimizer
fstar = fhandle(xstar); % true minimized value of f
kgk2 = 1e15; % length of the gradient vector
kpk2 = 1e15; % length of the vector p
alpha = 1.0; % step size for updating the solution, keeping it same for all iterations
k = 0; % iteration counter
while (k < maxits & kgk2 > gtol & kpk2 > xtol)
 [f,g,h] = fhandle(x); % function, gradient and Hessian
 kgk2 = sqrt(g(1)^2 + g(2)^2);
 p = -g/h; % solution of h*p = -g
 kpk2 = sqrt(p(1)^2 + p(2)^2);
 x = x + alpha*p; % update solution
 e = abs(f - fstar); % error
 fprintf("%d\t%4.2f\t%4.2f\n",k,f,e)
 %fprintf("%d\n",k)
 k = k + 1;
end

end
