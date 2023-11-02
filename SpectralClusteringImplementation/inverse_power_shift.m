function [lambda, y, iter] = inverse_power_shift(A, mu, tol, nmax, x0)

n = size(A, 1);
if nargin == 2
    tol = 1e-4;
    nmax = 10;
    x0 = rand(n, 1);
elseif nargin == 4
    x0 = rand(n, 1);
end

% initialization
lambda_old = 0;
iter = 0;

y = x0 / norm(x0);

[L, U] = lu(A - mu*eye(n));
z = L \ y; 
x = U \ z;
lambda = y'*x;
y = x / norm(x);
while iter < nmax && abs(lambda - lambda_old) >= tol * lambda
    lambda_old = lambda;
    z = L \ y;
    x = U \ z;
    lambda = y'*x;
    y = x / norm(x);
    iter = iter + 1;
end
fprintf("iter %i, diff %d \n", iter, abs(lambda - lambda_old))

lambda = 1 / lambda + mu;
