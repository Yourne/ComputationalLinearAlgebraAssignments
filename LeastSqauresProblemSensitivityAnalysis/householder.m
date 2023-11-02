function [W, R] = householder(A)

% R : upper triangular matrix
% W : householder reflector

[M, N] = size(A);
R = zeros(N, N);
W = zeros(M, N);

for k = 1 : N
    x = A(k:end, k);
    x(1) = x(1) + sign(x(1)) * norm(x); % equivale a v = sign(x(1))*e1 + x 
    x = x / norm(x); % normalized Householder reflector
    W(k:end, k) = x; 
    A(k:end, k:end) = A(k:end, k:end) - 2* x * (x' * A(k:end, k:end));
    R(1:k, k) = A(1:k, k);
end
end