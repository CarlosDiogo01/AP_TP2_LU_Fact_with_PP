function [A,L,U,P] = BLAS3LUPP(A,b)
% Block LU factorization with partial pivoting, overwriting L and U on A
% see ALGORITHM 2.10 in Applied Numerical Linear Algebra, J. Demmel, SIAM
% (2007), p.74. Size of blocks is b.
% Uses BLAS2LU to produce the LU dec. of rectangular submatrices of A

% Author Carlos Sá e Bruno Barbosa
% VERSION WITH PARTIAL PIVOTING

A_ORIGINAL = A;
n=length(A); P = eye(n);
[m n] = size(A);

start_time = tic;
for i=1:b:n-1
    
    for ln=i:i+b-1
        [~,p] = max(abs(A(ln:n,ln)));
    
        p = p+ln-1;
        if p~=ln
            A([ln p], :) = A([p ln], :);
            P([ln p], :) = P([p ln], :);
        end
    end
    %A(i+1:m,i)=A(i+1:m,i)/A(i,i);
    
    last=min(i+b-1,n);
    [A_res,~,~,P2]=BLAS2LUPP(A(i:m,i:last)); % step 1 (L22, L32)
    A(i:m,1:n) = P2 * A(i:m,1:n);
    P(i:m,1:n) = P2 * P(i:m,1:n);
    A(i:m,i:last) = A_res;
    
    % SIZE OF REMAINING BLOCK LARGER THAN b
    if n-i+1 > b  
        L22=tril(A(i:last,i:last),-1)+eye(b);
        A(i:last,i+b:n)=inv(L22)*A(i:last,i+b:n); % step 2 (U22)
        A(i+b:n,i+b:n)=A(i+b:n,i+b:n)-A(i+b:n,i:last)*A(i:last,i+b:n); % step 3 (U33)
    end
end
total_time = toc(start_time)

% Decomposition L U and P
L = tril(A,-1)+eye(size(A));
U = triu(A);

%Computing lu for validation
%LU_RESULT = lu(A_ORIGINAL)

% Compute Error
Relative_Error = norm(P*A_ORIGINAL - L*U)/norm(A)