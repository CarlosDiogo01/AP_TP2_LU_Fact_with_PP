function [ L,U,P ] = my_lu_piv( A )
%MY_LU_PIV Summary of this function goes here
%   Detailed explanation goes here
n = size(A,1); 
I = eye(n); 
O = zeros(n); 
L = I; 
U = O; 
P = I;

function change_rows(k,p)
    x = P(k,:); 
    P(k,:) = P(p,:); 
    P(p,:) = x;
    x = A(k,:); 
    A(k,:) = A(p,:); 
    A(p,:) = x;
    x = v(k); 
    v(k) = v(p); 
    v(p) = x;
end

function change_L(k,p)
    x = L(k,1:k-1); 
    L(k,1:k-1) = L(p,1:k-1);
    L(p,1:k-1) = x;
end

for k = 1:n
 if k == 1
     v(k:n) = A(k:n,k);
 else
    z = L(1:k-1,1:k -1)\ A(1:k-1,k);
    U(1:k-1,k) = z;
    v(k:n) = A(k:n,k)-L(k:n,1:k-1)*z;
 end
 
 if k<n
    x = v(k:n); p = (k-1)+find(abs(x) == max(abs(x))); % find index p
    change_rows(k,p);
    L(k+1:n,k) = v(k+1:n)/v(k);
    if k > 1
        change_L(k,p); 
    end
 end
 
 U(k,k) = v(k);

end

end


