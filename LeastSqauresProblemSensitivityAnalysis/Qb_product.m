function b = Qb_product(A, b, W)
% Implicit Q'*b product with reducted Q 
% A : coefficient matrix 
% b : constant vector
% W : householder reflector matrix
N = size(A, 2);
for k = 1 : N
    b(k:end) = b(k:end) - 2 * W(k:end, k) * (W(k:end, k)' * b(k:end));
end

end