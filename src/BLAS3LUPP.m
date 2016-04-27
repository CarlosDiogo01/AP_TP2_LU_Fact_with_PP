function [A,L,U,P] = BLAS3LUPP(A,b)
% Block LU factorization with partial pivoting, overwriting L and U on A
% see ALGORITHM 2.10 in Applied Numerical Linear Algebra, J. Demmel, SIAM
% (2007), p.74. Size of blocks is b.
% Uses BLAS2LU to produce the LU dec. of rectangular submatrices of A

% Author Carlos Sá e Bruno Barbosa
% VERSION WITH PARTIAL PIVOTING

A_ORIGINAL = A;
start_time = tic;
n=length(A);
P = eye(n);

for i=1:b:n-1
    % apply row permutations to A and L
    last=min(i+b-1,n);
    % locate pivot's position lu(A)
    [~,p] = max(abs(A(i:last,i)));
    p = p+i-1;
    % swap rows only if pivot is not i
    if p~=i
        A([i p], :) = A([p i], :);
        P([i p], :) = P([p i], :);
    end
    % BLAS2LUPP call
    [A(i:n,i:last),~,~,P2]=BLAS2LUPP(A(i:n,i:last)); % step 1 (L22, L32)
    % swap lines
    A(i:last,1:i-1) = P2*A(i:last,1:i-1); % left side of the block
    A(i:last,i+b:n) = P2*A(i:last,i+b:n); % right side of the block
    P(i:last,:) = P2*P(i:last,:); % updates local permutations matrix
    
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
LU_RESULT = lu(A_ORIGINAL)
% Compute Error
Relative_Error = norm(P*A_ORIGINAL - L*U)/norm(A)
