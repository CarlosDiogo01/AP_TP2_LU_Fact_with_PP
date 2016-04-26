function [A L U] = BLAS2LU(A)
% LU factorization with partial pivoting, overwriting L and U on A
% see ALGORITHM 2.9 in Applied Numerical Linear Algebra, J. Demmel, SIAM
% (2007), p.72
% This differs from LUfact1 because A may not be square, has m rows and n
% columns, with m>=n
%
% VERSION WITHOUT PARTIAL PIVOTING

A_ORIGINAL = A;

start_time = tic;
[m n]=size(A);
for i=1:min(m-1,n)
    % apply row permutations to A and L
    A(i+1:m,i)=A(i+1:m,i)/A(i,i);
    if i<n
        A(i+1:m,i+1:n)=A(i+1:m,i+1:n)-A(i+1:m,i)*A(i,i+1:n);
    end
end
total_time = toc(start_time)

% Decomposition L U and P
L = tril(A); 
for j=1:n
    L(j,j) = 1;
end
U = triu(A);

% Compute Error
Relative_Error = norm(A_ORIGINAL - L*U)/norm(A)
