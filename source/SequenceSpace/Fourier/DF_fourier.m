function mat = DF_fourier(a,order, params)
% when computing DF_k/Da_j, we differentiate with respect to a_j and a_{-j}
order = max(size(a))-1;
mat=zeros(order+1, order+1);
a=[a, zeros(1,size(a,2))];
duba=dubconvvec(a,a, order);

L = params.L*2;
mu = params.mu;
nu = params.nu;

for k = 0:order 
    mat(k+1,1) = -3*duba(k+1)  + 2*nu*a(k+1);
    for j = 1:k-1
        mat(k+1,j+1) = -3*duba(k-j+1)  + 2*nu*a(k-j+1);
        if abs(k+j)<order+1
            mat(k+1,j+1)=mat(k+1,j+1)-3*duba(k+j+1)+2*nu*a(k+j+1);
        end
    end
    
    mat(k+1,k+1) = -mu-(1-(k*pi/L)^2)^2 - 3*duba(1)+ 2*nu*a(1);
    if (2*k < order+1) && (k > 0)
        mat(k+1,k+1) = mat(k+1,k+1) - 3*duba(2*k+1) + 2*nu*a(2*k+1);
    end
    
    for j = k+1:order
        mat(k+1,j+1) = -3*duba(abs(k-j)+1) + 2*nu*a(abs(k-j)+1);
        if abs(k+j)<order+1
            mat(k+1,j+1)=mat(k+1,j+1)-3*duba(k+j+1)+2*nu*a(k+j+1);
        end
    end
end
end
