function [Q, R] = mgs(A)
% reduced orthogonalization with modified gram schmidt method

[~, N] = size(A);
Q = A;
R = zeros(N, N);

for i = 1:N
    R(i, i) = norm(Q(:, i));
    Q(:, i) = Q(:, i) / R(i, i);
    for j = i+1:N
        R(i, j) = Q(:,i)' * Q(:, j);
        Q(:, j) = Q(:,j) - R(i, j)*Q(:,i);
    end
end

