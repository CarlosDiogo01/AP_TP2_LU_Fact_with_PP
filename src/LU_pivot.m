function [L,U,P]=LU_pivot(A)
% LU factorizatin wiht partial pivoting
[m n]=size(A);
L=eye(n); P=L; U=A;
for i=1:min(m-1,n)
    [~,p]=max(abs(U(i:n,i)));
    p=p+i-1;
    if p~=i
        % swap rows p and i in U
        U([i p], :) = U([p i], :);
        % swap rows p and i in P
        P([i p], :) = P([p i], :);
        if i >= 2
             L([i p],1:i-1) = L([p i],1:i-1);
        end
    end
    L(i+1:m,i) = U(i+1:m,i) / U(i,i);
    U(i+1:m,:) = U(i+1:m,:) - L(i+1:m,i)*U(i,:);
end