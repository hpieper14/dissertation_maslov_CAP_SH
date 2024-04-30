function val=mfld_poly(params,coeff, Q,Lambda)
    order=params.mfld.order;
    % Calculate the coefficient K 
    maxKmn=((order+1)*abs(real(Lambda(1))) - abs(Lambda(1)))^(-1);
 
    K_N = max(abs(Q),[],'all')*max(abs(Q^(-1)),[],'all')*maxKmn;
    
    % extract the coefficients a_{mn}
    
    a=coeff(1:order+1,1:order+1,1);
    
    A=zeros(order+1, 2*order);
    B=zeros(2*order, order+1);
    C=zeros(2*order, 2*order);
    
    bara=[a, A ; B, C];
    %bara=coeff(1:3*order+1, 1:3*order+1,1);
    
    disp('Calculating Y0.')

    %%%%%%
    % Y0 %
    %%%%%%

    quadsum=0;
    suborder=order+1;
    while suborder < 2*order+1
        for i=0:suborder
            j=suborder-i;
            astarij=starhat(bara, bara, i,j);
            quadsum=quadsum+abs(astarij);
        end
        suborder=suborder+1;
    end
   
    cubsum=0;
    suborder=order+1;
    while suborder<3*order+1
        for i=0:suborder
            j=suborder-i;
            tripstara=tripstarhat(bara,bara,bara,i,j);
            cubsum=cubsum+abs(tripstara); 
        end
        suborder=suborder+1;
    end
    
    Y0 = K_N*(params.nu*quadsum+cubsum);
    disp(Y0)
    disp('Calculating Z1.')
    
    %%%%%%
    % Z1 %
    %%%%%%
    
    lowerquadsum=0;
    suborder=1;
    while suborder < order+1
        for i=0:suborder
            j=suborder-i;
            astarij=starhat(bara, bara, i,j);
            lowerquadsum=lowerquadsum+abs(astarij);
        end
        suborder=suborder+1;
    end
    
    linsum=0;
    suborder=1;
    while suborder<order+1
        for i=0:suborder
            j=suborder-i;
            linsum=linsum+abs(bara(i+1,j+1));
        end
        suborder=suborder+1;
    end
    
    Z1= K_N*(2*params.nu*linsum + 3*lowerquadsum);
    disp(Z1);
    disp('Calculating Z2.')
    %%%%%%
    % Z2 %
    %%%%%%
    
    Z2=@(r)K_N*(6*linsum+2*params.nu+3*r);
    disp(Z2)
    
    b=K_N*3;
    a=K_N*(6*linsum + 2*params.nu);
  
    disp(a)
    disp(b)
    % this section of code builds an interval on which the sup of the polynomial is negative 
    poly=@(r)(Z1+Z2(r)*r)*r+Y0-r;
    
    %figure(3)
    %fplot(@(r) poly(r));
    
    p=[b, a, Z1-1, Y0];
    
    val=roots(p);
end
