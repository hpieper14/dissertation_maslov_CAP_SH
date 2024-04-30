% x = (phi1, phi2, psi, a1, a2, a3, a4)
% store ai as row vectors
function fun = FHomoclinic(x, mflds, params)

    psi1 = params.rho*cos(x.psi);
    psi2 = params.rho*sin(x.psi);
    
    P = real(mfld_one_point(psi1,psi2,mflds.stable.coeffs, params));
    Q = real(mfld_one_point(x.phi1, x.phi2, mflds.unstable.coeffs, params));
   
    
    fun = cell(7,1);
    
    %%%%%%%%%%%%%%%
    % f_1 to f_3  %
    %%%%%%%%%%%%%%%
    
    fun{1} = x.a1(1) + 2*sum(x.a1(2:end)) - P(1);
    fun{2} = x.a2(1) + 2*sum(x.a2(2:end)) - P(2);
    fun{3} = x.a4(1) + 2*sum(x.a4(2:end)) - P(4); 
    
    %%%%%%%%%%%%%%%
    % f_4 to f_6  %
    %%%%%%%%%%%%%%%
    
    m = params.cheb.order;
    pm_ones = (-ones(1,m)).^(1:m);
    
    if max(size(x.a1)) == m
        x.a1=[x.a1, 0];
        x.a2=[x.a2, 0];
        x.a3=[x.a3, 0];
        x.a4=[x.a4, 0];
    end
    
    
    f4to6=zeros(3,m);
  
    f4to6(1,1) = x.a1(1) + 2*sum(pm_ones.*x.a1(2:end)) - Q(1);
    f4to6(2,1) = x.a2(1) + 2*sum(pm_ones.*x.a2(2:end)) - Q(2);
    f4to6(3,1) = x.a3(1) + 2*sum(pm_ones.*x.a3(2:end)) - Q(3);
    
    for i = 1:m-1
        f4to6(1,i+1) = 2*i*x.a1(i+1) - params.L*(x.a2(i)-x.a2(i+2));
        f4to6(2,i+1) = 2*i*x.a2(i+1) - params.L*(x.a3(i)-x.a3(i+2));
        f4to6(3,i+1) = 2*i*x.a3(i+1) - params.L*(x.a4(i)-x.a4(i+2));
    end
    
    fun{4} = f4to6(1,:);
    fun{5} = f4to6(2,:);
    fun{6} = f4to6(3,:);
    
    %%%%%%%%
    % f_7  %
    %%%%%%%%
    
    f7=zeros(1,m);
    f7(1) = x.a4(1) + 2*sum(pm_ones.*x.a4(2:end)) - Q(4);
    
    c4=zeros(1,m+1);
    a1a1=chebstar2(x.a1,x.a1,m+1);
    a1a1a1=chebstar3(x.a1,x.a1,x.a1,m+1);
    
    for i=1:m+1
        c4(i) = -2*x.a3(i) - (1+params.mu)*x.a1(i) + params.nu*a1a1(i) ...
            - a1a1a1(i);
    end
    
    for i = 1:m-1
        f7(i+1) = 2*i*x.a4(i+1) - params.L*(c4(i) - c4(i+2));
    end
    
    fun{7} = f7;
  
end
