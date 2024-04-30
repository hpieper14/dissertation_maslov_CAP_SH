function coeff = calc_proj_coeff(eigenvalues,eigenvectors,params)
    order=params.mfld.order;
    coeff=zeros(order+1,order+1,4);

    

    e1=eigenvectors(:,1)*params.scale;
    e2=eigenvectors(:,2)*params.scale;
    lam1=eigenvalues(1);
    lam2=eigenvalues(2);

    disp(order)
    
    Df0=JacSH(0,params.mu,params.nu); 
    % the zeroth order coefficient is the equilibrium, corresponding to 0.
    coeff(2,1,:)=e1;
    coeff(1,2,:)=e2;
    
    % suborder=m+n and corresponds to the (m,n)th coefficient
   suborder=2;
   while suborder < order + 1 
        for j=suborder:-1:0
            i=suborder-j;
            Aij=(Df0-(i*lam1+j*lam2)*eye(4))^(-1);
            a=coeff(1:i+1,1:j+1,1);
            staraaij=starhat(a,a,i,j);
            staraaaij=tripstarhat(a,a,a,i,j);
            B=[0;0;0;-params.nu*staraaij+staraaaij];
            pij=Aij*B;
            coeff(i+1,j+1,:)=pij;
        end
        suborder=suborder+1;
   end
end