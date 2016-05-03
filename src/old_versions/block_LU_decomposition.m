function [L,U] = block_LU_decomposition(A)

n_group = size(A,1);
for k = 1 : n_group
    %[L_tmp, U_tmp] = lu(A{k,k});
    [L_tmp, U_tmp] = lu(A);
    L{k,k} = L_tmp;
    U{k,k} = U_tmp;
    
    for i = k+1 : n_group
        L(i,k) = A(i,k)/U_tmp;
        U(k,i) = L_tmp\A(k,i);
    end  
    for j=k+1 : n_group
        for i=k+1 : n_group
            A{i,j} = A{i,j} - L{i,k}*U{k,j};
        end
    end
end