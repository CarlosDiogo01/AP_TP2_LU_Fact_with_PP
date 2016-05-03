function [A,L,U] = BLAS3LU(A,b)
% Block LU factorization with partial pivoting, overwriting L and U on A
% see ALGORITHM 2.10 in Applied Numerical Linear Algebra, J. Demmel, SIAM
% (2007), p.74. Size of blocks is b.
% Uses BLAS2LU to produce the LU dec. of rectangular submatrices of A

% VERSION WITHOUT PARTIAL PIVOTING %

A_ORIGINAL = A;

start_time = tic;
n=length(A);
for i=1:b:n-1
    % apply row permutations to A and L
    last=min(i+b-1,n);
    % BLAS2LUPP call
    A(i:n,i:last)=BLAS2LU(A(i:n,i:last)); % step 1 (L22, L32)
    
    
    
    
    if n-i+1 > b  % SIZE OF REMAINING BLOCK LARGER THAN b
        L22=tril(A(i:last,i:last),-1)+eye(b);
        A(i:last,i+b:n)=inv(L22)*A(i:last,i+b:n); % step 2 (U22)
        A(i+b:n,i+b:n)=A(i+b:n,i+b:n)-A(i+b:n,i:last)*A(i:last,i+b:n); % step 3 (U33)
    end
end
total_time = toc(start_time)

% Decomposition L U and P
L = tril(A); 
for j=1:n
    L(j,j) = 1;
end
U = triu(A);
%Computing lu for validation
% LU_RESULT = lu(A_ORIGINAL)

% Compute Error
Relative_Error = norm(A_ORIGINAL - L*U)/norm(A)
