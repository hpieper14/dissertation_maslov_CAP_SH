% Jacobian matrix for fourier series of ODE with Neumann BC
function mat = DFFourierNBC(a,order, params,side)
% when computing DF_k/Da_j, we differentiate with respect to a_j and a_{-j}
order = max(size(a))-1;
mat=zeros(order+1, order+1);

a=[a, zeros(1,size(a,2))];
duba=dubconvvec(a,a, 2*order,side);

L = params.L;
mu = params.mu;
nu = params.nu;

% kth row corresponds to F_k 
% jth column corresponds to d/da_j
for k = 0:order 
    mat(k+1,1) = -3*duba(k+1)  + 2*nu*a(k+1);
    if k == 0
        mat(k+1,1) = mat(k+1,1) - mu - 1;
    end
    for j = 1:order       
        % df_k/daj with k == j
        if k == j
            mat(k+1,j+1) = -mu-(1-(k*pi/L)^2)^2 - 3*(duba(abs(k-j)+1)+duba(abs(k+j)+1))...
                + 2*nu*(a(abs(k-j)+1)+a(abs(k+j)+1));
        % df_k/daj with k =\= j
        else
            mat(k+1, j+1) = -3*(duba(abs(k-j)+1)+duba(abs(k+j)+1)) ...
                + 2*nu*(a(abs(k-j)+1) + a(abs(k+j)+1));
        end
     end
end
end
