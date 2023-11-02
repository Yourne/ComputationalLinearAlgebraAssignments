function lambda = power_iter_deflate(A, k)

% A : N x N symmetric real diagonally dominant matrix
% k : number of eigenvalues to compute
% lambda : top k largest eigevalues of A

n = size(A, 1);
if n-k < 2
    error("this function computes up to " + n-2 + " eigenvalues");
end

lambda = zeros(k, 1);

[lambda(1), w] = power_iteration(A);

for j = 2:k
    [lambda(j), w, A] = power_deflation(A, lambda(j-1), w);
end

end