function [A,L,U,P] = BLAS2LUPP(A)
% LU factorization with partial pivoting, overwriting L and U on A
% see ALGORITHM 2.9 in Applied Numerical Linear Algebra, J. Demmel, SIAM
% (2007), p.72
% This differs from LUfact1 because A may not be square, has m rows and n
% columns, with m>=n

% Author Carlos Sá e Bruno Barbosa
% VERSION WITH PARTIAL PIVOTING

A_ORIGINAL = A;
start_time = tic;
[m n]=size(A);
P = eye(n);
for i=1:min(m-1,n)
%     A(i:m,i) = A(i:m,i) - A(i:m,1:i-1) * A(1:i-1,i);
    % apply row permutations to A and L
    % locate pivot's position lu(A)
    [~,p] = max(abs(A(i:n,i)));
    p = p+i-1;
    % swap rows only if pivot is not i
    if p~=i
        A([i p], :) = A([p i], :);
        P([i p], :) = P([p i], :);
    end
    A(i+1:m,i)=A(i+1:m,i)/A(i,i);
    if i<n
        A(i+1:m,i+1:n)=A(i+1:m,i+1:n)-A(i+1:m,i)*A(i,i+1:n);
    end
end

%Measuring Time
total_time = toc(start_time);

% Decomposition L U and P
L = tril(A); 
for j=1:n
    L(j,j) = 1;
end
U = triu(A);

%Computing lu for validation
%LU_RESULT = lu(A_ORIGINAL);

% Compute Error
% Relative_Error = norm(P*A_ORIGINAL - L*U)/norm(A)
