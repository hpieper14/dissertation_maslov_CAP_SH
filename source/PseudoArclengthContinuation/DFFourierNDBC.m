function mat = DFFourierNDBC(a,order, params)
% when computing DF_k/Da_j, we differentiate with respect to a_j and a_{-j}
mat=zeros(2*order+1, 2*order+1);

if (size(a,1) == 1) == 0
    a = a';
end

side = 2;

a=[zeros(1,order), a, zeros(1,order)];
duba=dubconvvec(a,a, 2*order, side);

L = params.fourier.L;
mu = params.mu;
nu = params.nu;

% lth row corresponds to F_l
% kth column corresponds to d/da_k


for l = -order:order 
    matrowLind = l  + order + 1; 
    for k = -order:order
        matcolKind = k + order + 1; 
        aLKind = l - k + 2*order+1;
        mat(matrowLind, matcolKind) = -3*(duba(aLKind)) ...
                + 2*nu*(a(aLKind));
        if k == l
            mat(matrowLind,matcolKind) = mat(matrowLind,matcolKind)-mu-(1-(k*pi/L)^2)^2;
        end
    end
end
        
