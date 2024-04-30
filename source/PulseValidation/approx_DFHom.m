function approx_df = approx_DFHom(x, epsilon,params,mflds)
    F=FHomoclinic(x,mflds, params);
    m=params.cheb.order;
    approx_df=zeros(4*m+3);
    
    totalF=[F{1}, F{2}, F{3}, F{4}, F{5}, F{6}, F{7}];
    
    for i=1:4*m+3
        perturb=zeros(1,4*m+3);
        perturb(i) = epsilon;
        totalx=[x.phi1,x.phi2,x.psi,x.a1,x.a2,x.a3,x.a4];
        perturbx = totalx+perturb;
        newx.phi1=perturbx(1);
        newx.phi2=perturbx(2);
        newx.psi=perturbx(3);
        newx.a1=perturbx(4:m+3);
        newx.a2=perturbx(3+m+1:2*m+3);
        newx.a3=perturbx(3+2*m+1:3*m+3);
        newx.a4=perturbx(3+3*m+1:4*m+3);
        
        perturbF=FHomoclinic(newx,mflds,params);
        
        totalPerturb=[perturbF{1}, perturbF{2}, perturbF{3}, perturbF{4}, ...
            perturbF{5}, perturbF{6}, perturbF{7}];
        diff=totalPerturb-totalF;
        approx_df(:,i)=(diff')./epsilon;
    end
end
