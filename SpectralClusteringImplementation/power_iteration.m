function [lambda, x] = power_iteration(A, tol, nmax, x0)

% compute the largest eigenvalue of the matrix in argument
% A matrix must be real and nonsingular
% x is an arbitrary initial vector column vector

N = size(A);
if N(1) ~= N(2)
    error ('Only for square matrices'); 
end

if nargin == 1
    x0 = rand(N(1), 1);
    tol = 1e-8;
    nmax = 1000;
elseif nargin == 3
    x0 = rand(N(1), 1);
end


x = x0 / norm(x0);
x = A * x;
lambda_old = 0;
lambda =  x.' * A * x;
n = 0;
while (n < nmax && abs(lambda - lambda_old) >= tol*abs(lambda))
    lambda_old = lambda;
    x = A * x;
    x = x / norm(x);
    lambda = x.' * A * x;
    n = n + 1;
end
fprintf("n iter: %d \n", n)
end