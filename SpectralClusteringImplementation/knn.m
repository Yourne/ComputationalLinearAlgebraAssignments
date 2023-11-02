function W = knn(X, k)

% X : m by n data matrix
% k : nearest neighbor parameter
% W : nearest neighbor matrix

% compute the similary matrix 
% similary function: gaussian kernel
[n, m] = size(X);
X = reshape(X, [n, 1, m]);
S = (X - pagetranspose(X)).^2;
S = sqrt(sum(S, 3));

% compute the adjacency matrix

s_row = sort(S, 2); 
s_col = sort(S, 1);

% qui andrà scelta una tolleranza eps
% idx_row = (abs(S - s_row(:, k+1))) <= eps;
% idx_col = (abs(S - s_col(k+1, :))) <= eps;
idx_row = S <= s_row(:, k+1);
idx_col = S <= s_col(k+1, :);
idx = or(idx_row, idx_col); 
W = zeros(n, n);
W(idx) = S(idx);

end