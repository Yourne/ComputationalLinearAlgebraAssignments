% build hilbert matrix
M = 100; N = 6;
H = hilb(M);
H = H(:, 1:N);

x_ = [1:N]';
b = H * x_;

% estimate x with MATLAB solver
x = H \ b;

% least squares conditioning parameters
theta = acos(norm(H*x) / norm(b));
eta = norm(H) * norm(x) / norm(H*x);
cond_H2x = (cond(H) + cond(H)^2 * tan(theta) / eta)*2^-53;
fprintf("QR accuracy should be less than %e \n", cond_H2x)

% Householder solution
[W, R] = householder(H);
y = Qb_product(H, b, W);
x = R \ y(1:N);
fprintf("Householder accuracy  %e \n", norm(x - x_) / norm(x_))

% MGS solution
[Q, R] = mgs(H); 
x = R \ (Q' * b);
fprintf("MGS accuracy  %e \n", norm(x - x_) / norm(x_))
fprintf("orthogonality test Q %e \n", norm(Q'*Q - eye(N)))
