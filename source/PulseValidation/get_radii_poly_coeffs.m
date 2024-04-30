function [Y,Z, Z0] = get_radii_poly_coeffs(nu,x,mflds,params)

    m = params.cheb.order;
    N = params.mfld.order;
    
    norms = get_Aa_norms(nu, x, params, mflds);

    disp('Computing Y.')

    F = FHomoclinic(x, mflds, params);
    F = [F{1}; F{2}; F{3}; F{4}'; F{5}'; F{6}'; F{7}'];
    
    delta_s = mflds.stable.error;
    delta_u = mflds.unstable.error;

    F(1:7) = F(1:7)+[delta_s*ones(3,1);delta_u*ones(4,1)];
    
    DF = DF_homoclinic(x, params, mflds);
    Am = DF^(-1);

    vF_1=abs(Am(1,:)*F);
    vF_2=abs(Am(2,:)*F);
    vF_3=abs(Am(3,:)*F);
    vF_4=abs(Am(3+1:3+m,:)*F);
    vF_5=abs(Am(3+m+1:3+2*m,:)*F);
    vF_6=abs(Am(3+2*m+1:3+3*m,:)*F);
    vF_7=abs(Am(3+3*m+1:3+4*m,:)*F);

    Y=zeros(7,1);
    
    nu_power=nu.^(0:m-1);

    Y(1)=vF_1;
    Y(2)=vF_2;
    Y(3)=vF_3;
    Y(4)=sum(vF_4.*nu_power');
    Y(5)=sum(vF_5.*nu_power');
    Y(6)=sum(vF_6.*nu_power');
    Y(7)=sum(vF_7.*nu_power');

    disp('The bound Y is: ')
    disp(Y)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Computation of the bound Z %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%
    %%%% Z0 %%%%
    %%%%%%%%%%%%
    
    disp('Computing Z0.')

    % We compute the norms of the sub operators of the operator B
    B=eye(4*m+3)-Am*DF;
    B_norms=zeros(7,7);
    B_norms(1:3,1:3)=abs(B(1:3,1:3));
    
    for j=4:7
        B_norms(1:3,j)=max(abs(B(1:3,4*(0:m-1)+j).*repmat(nu.^-(0:m-1),3,1)),[],2);
    end
    
    for i=4:7
        B_norms(i,1:3)=sum(abs(B(4*(0:m-1)+i,1:3)).*repmat(nu_power',1,3),1);
    end
    
    K_B=zeros(7,7);
    for i=4:7
        for j=4:7
            K_B(i,j)=max(mag(sum(abs(B(4*(0:m-1)+i,4*(0:m-1)+j)).*repmat(nu_power,m,1)').*(nu.^-(0:m-1))));
        end
    end
    B_norms(4:7,4:7)=max(K_B(4:7,4:7),0);

    % We compute the bound Z0
    Z0=zeros(7,1);
    Z0(1:7)=sum(B_norms(1:7,1:7),2);

    disp('The bound Z0 (linear in r) is ')
    disp(Z0)

    %%%%%%%%%%%%%%%%
    %%%% LAMBDA %%%%
    %%%%%%%%%%%%%%%%

    disp('Computing Lambda (manifolds) and verifying valid domain.')

    % We choose ru_star and rs_star to be small enough.
    
    
    %%% THIS SECTION NEEDS ATTENTION -- this is ensuring that we say in te domain of definition in the unstable manifold
    %%% WHY CAN WE USE phi1 and phi2 instead of phi1 - iphi2 and phi1 +
    %%% iphi2???
    
    ru_star=1e-4;
    rs_star=1e-4;

    sigma_s=infsup(x.psi-rs_star,x.psi+rs_star);
    Psi_1=params.rho*cos(sigma_s);
    Psi_2=params.rho*sin(sigma_s);
    
    a1 = x.phi1 + 1i*x.phi2;
    a2 = x.phi1 - 1i*x.phi2;
    
    sigma_u1 = infsup(x.phi1-ru_star, x.phi1+ru_star);
    sigma_u2 = infsup(x.phi2-ru_star, x.phi2+ru_star);

    hat_sigma_u1=intval(mag(sigma_u1));
    hat_sigma_u2=intval(mag(sigma_u2));
    
    disp('Now verifying that we have not left the domain of definition of the unstable manifold.');
    disp('')
    disp(['||hat_sigma_u1|| = ',num2str((mag(hat_sigma_u1))),' ||hat_sigma_u2|| = ',num2str((mag(hat_sigma_u2))) ])

    if sup(hat_sigma_u1)>=1 || sup(hat_sigma_u2)>=1
        disp('Stop! The parameter values are not within the defined region of the unstable manifold.')
        return
    else 
        disp('Good to go!')
    end


    % We compute the bounds LAMBDA^{(s,i)}
    
    
    powers_s11 = zeros(N+1, N+1);
    powers_s21 = zeros(N+1, N+1);
    powers_s22 = zeros(N+1, N+1);
    powers_s1 = zeros(N+1, N+1);
    powers_s2 = zeros(N+1, N+1);
    
    
    powers_u11 = zeros(N+1, N+1);
    powers_u21 = zeros(N+1, N+1);
    powers_u22 = zeros(N+1, N+1);
    
    a1 = params.rho*cos(x.psi) + 1i*params.rho*sin(x.psi);
    a2 = params.rho*cos(x.psi) - 1i*params.rho*sin(x.psi); 
    
    b1 = x.phi1 + 1i*x.phi2;
    b2 = x.phi1 - 1i*x.phi2;
    
    %a1 = Psi_1 + 1i*Psi_2;
    %a2 = Psi_1 - 1i*Psi_2; 
    
    %b1 = sigma_u1;
    %b2 = sigma_u2;
    
    for q = 0:N 
        for n = 0:N 
            powers_s11(q+1,n+1) = q*(q-1)*a1^(q-1)*a2^n + 2*q*n*a1^(q-1)*a2^(n-1)...
              +n*(n-1)*a1^q*a2^(n-2);
            powers_s21(q+1,n+1) = 1i*q*(q-1)*a1^(q-2)*a2^n - 1i*n*(n-1)*a1^q*a2^(q-2);
            powers_s22(q+1, n+1) = -q*(q-1)*a1^(q-2)*a2^n+2*q*n*a1^(q-1)*a2^(n-1)...
                -n*(n-1)*a1^q*a2^(n-2);
            
            powers_s1(q+1, n+1) = q*a1^(q-1)*a2^n + n*a1^q*a2^(n-1);
            powers_s2(q+1, n+1) = 1i*q*a1^(q-1)*a2^n - 1i*n*a1^q*a2^(n-1);
            
            powers_u11(q+1,n+1) = q*(q-1)*b1^(q-1)*b2^n + 2*q*n*b1^(q-1)*b2^(n-1)...
              +n*(n-1)*b1^q*b2^(n-2);
            powers_u21(q+1,n+1) = 1i*q*(q-1)*b1^(q-2)*b2^n - 1i*n*(n-1)*b1^q*b2^(n-2);
            powers_u22(q+1, n+1) = -q*(q-1)*b1^(q-2)*b2^n+2*q*n*b1^(q-1)*b2^(n-1)...
                -n*(n-1)*b1^q*b2^(n-2);
            
        end
    end
    
    p1 = mflds.stable.coeffs(:,:,1);
    p2 = mflds.stable.coeffs(:,:,2);
    p4 = mflds.stable.coeffs(:,:,4);
    
    SUM_S = zeros(4,5);
    
    SUM_S(1,1)=sum(sum(powers_s11.*abs(p1)));
    SUM_S(2,1)=sum(sum(powers_s11.*abs(p2)));
    SUM_S(4,1)=sum(sum(powers_s11.*abs(p4)));

    
    SUM_S(1,2)=sum(sum(powers_s21.*abs(p1)));
    SUM_S(2,2)=sum(sum(powers_s21.*abs(p2)));
    SUM_S(4,2)=sum(sum(powers_s21.*abs(p4)));

    
    SUM_S(1,3)=sum(sum(powers_s22.*abs(p1)));
    SUM_S(2,3)=sum(sum(powers_s22.*abs(p2)));
    SUM_S(4,3)=sum(sum(powers_s22.*abs(p4)));

    SUM_S(1,4)=sum(sum(powers_s1.*abs(p1)));
    SUM_S(2,4)=sum(sum(powers_s1.*abs(p2)));
    SUM_S(4,4)=sum(sum(powers_s1.*abs(p4)));

    SUM_S(1,5)=sum(sum(powers_s2.*abs(p1)));
    SUM_S(2,5)=sum(sum(powers_s2.*abs(p2)));
    SUM_S(4,5)=sum(sum(powers_s2.*abs(p4)));

    LAMBDA_S_1=abs(SUM_S(1,1))*params.rho^2+2*abs(SUM_S(1,2))*params.rho^2 ...
        +abs(SUM_S(1,3))*params.rho^2+abs(SUM_S(1,4))*params.rho+abs(SUM_S(1,5))*params.rho;
    LAMBDA_S_2=abs(SUM_S(2,1))*params.rho^2+2*abs(SUM_S(2,2))*params.rho^2 ...
        +abs(SUM_S(2,3))*params.rho^2+abs(SUM_S(2,4))*params.rho+abs(SUM_S(2,5))*params.rho;
    LAMBDA_S_4=abs(SUM_S(4,1))*params.rho^2+2*abs(SUM_S(4,2))*params.rho^2 ...
        +abs(SUM_S(4,3))*params.rho^2+abs(SUM_S(4,4))*params.rho+abs(SUM_S(4,5))*params.rho;

    % We compute the bounds LAMBDA^{(u,i)}
    
    q1 = mflds.unstable.coeffs(:,:,1);
    q2 = mflds.unstable.coeffs(:,:,2);
    q3 = mflds.unstable.coeffs(:,:,3);
    q4 = mflds.unstable.coeffs(:,:,4);
    
    SUM_U = zeros(4,4);
    
    SUM_U(1,1)=sum(sum(powers_u11.*abs(q1)));
    SUM_U(2,1)=sum(sum(powers_u11.*abs(q2)));
    SUM_U(3,1)=sum(sum(powers_u11.*abs(q3)));
    SUM_U(4,1)=sum(sum(powers_u11.*abs(q4)));

    SUM_U(1,2)=sum(sum(powers_u21.*abs(q1)));
    SUM_U(2,2)=sum(sum(powers_u21.*abs(q2)));
    SUM_U(3,2)=sum(sum(powers_u21.*abs(q3)));
    SUM_U(4,2)=sum(sum(powers_u21.*abs(q4)));

    SUM_U(1,3)=sum(sum(powers_u22.*abs(q1)));
    SUM_U(2,3)=sum(sum(powers_u22.*abs(q2)));
    SUM_U(3,3)=sum(sum(powers_u22.*abs(q3)));
    SUM_U(4,3)=sum(sum(powers_u22.*abs(q4)));

    LAMBDA_U_1=abs(SUM_U(1,1))+2*abs(SUM_U(1,2))+abs(SUM_U(1,3));
    LAMBDA_U_2=abs(SUM_U(2,1))+2*abs(SUM_U(2,2))+abs(SUM_U(2,3));
    LAMBDA_U_3=abs(SUM_U(3,1))+2*abs(SUM_U(3,2))+abs(SUM_U(3,3));
    LAMBDA_U_4=abs(SUM_U(4,1))+2*abs(SUM_U(4,2))+abs(SUM_U(4,3));

    %%%%%%%%%%%%%%%
    %%%% Z_INF %%%%
    %%%%%%%%%%%%%%%

    disp('Computing z_inf.')

    z_inf=zeros(7,1);

    z_inf(4)=(params.L/(2*m))*(nu+1/nu);
    z_inf(5)=(params.L/(2*m))*(nu+1/nu);
    z_inf(6)=(params.L/(2*m))*(nu+1/nu);
    z_inf(7)=(params.L/(2*m))*(nu+1/nu)*(8*params.nu*norms.a1+48*norms.a1^2 +3 + params.mu);
    
    disp('The bound z_inf is ')
    disp(z_inf(4:7))

    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%% LAST PART OF Z %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%

    % We compute the FFT sums
 
    vI = [zeros(1,m),nu.^-(m:3*m)];
    x.a1 = [x.a1, zeros(1,2*m)];

    [a1vI]=chebstar2(abs(x.a1),vI, 2*m);
    [a1vIvI]=chebstar3(abs(x.a1),vI,vI, 3*m); 

    
    % We compute the last part of the bound Z
    disp('Computing the last part of Z.')

    cte_s=(4*pi/(log(1/params.rho)))*delta_s*params.rho+2/(nu^m);
    cte_u=(4*pi/(log(1/(1-ru_star))))*delta_u+2/(nu^m);

    Z=zeros(7,3);
    
     Z(1:3,1)=cte_s*sum(norms.A(1:3,1:3),2)+cte_u*sum(norms.A(1:3,4:7),2)+...
        +params.L*abs(Am(1:3, 3+3*m+1:3+4*m))*(2*params.nu*a1vI(1:m)'+...
        2*params.nu*a1vI(3:m+2)' + 3*a1vIvI(1:m)' +3*a1vIvI(3:m+2)');
    
    
        
      Z(4,1)=cte_s*sum(norms.A(4,1:3),2) ... 
        +cte_u*sum(abs(Am(4,4)) + abs(Am(4,3+m+1)) + abs(Am(4,3+2*m+1))+ abs(Am(4,3+3*m+1)))...
        +cte_u*sum(abs(Am(4:3+m,3+3*m+1)))...
        +params.L*abs(Am(4, 3+3*m+1:3+4*m))*((2*params.nu*a1vI(1:m)'+...
        2*params.nu*a1vI(3:m+2)' + 3*a1vIvI(1:m)' +3*a1vIvI(3:m+2)').*nu_power')...
        + z_inf(4);
    

    Z(5,1)=cte_s*sum(norms.A(5,1:3),2) ... 
        +cte_u*sum(abs(Am(3+m+1,4)) + abs(Am(3+m+1,3+m+1)) + abs(Am(3+m+1,3+2*m+1))+ abs(Am(3+m+1,3+3*m+1)))...
        +cte_u*sum(abs(Am(3+m+1:3+2*m,3+3*m+1)))...
        +params.L*abs(Am(3+m+1, 3+3*m+1:3+4*m))*((2*params.nu*a1vI(1:m)'+...
        2*params.nu*a1vI(3:m+2)' + 3*a1vIvI(1:m)' +3*a1vIvI(3:m+2)').*nu_power')...
        + z_inf(5);
    
     Z(6,1)=cte_s*sum(norms.A(6,1:3),2) ... 
        +cte_u*sum(abs(Am(3+2*m+1,4)) + abs(Am(3+2*m+1,3+m+1)) + abs(Am(3+2*m+1,3+2*m+1))+ abs(Am(3+2*m+1,3+3*m+1)))...
        +cte_u*sum(abs(Am(3+2*m+1:3+3*m,3+3*m+1)))...
        +params.L*abs(Am(3+2*m+1, 3+3*m+1:3+4*m))*((2*params.nu*a1vI(1:m)'+...
        2*params.nu*a1vI(3:m+2)' + 3*a1vIvI(1:m)' +3*a1vIvI(3:m+2)').*nu_power')...
        + z_inf(6);
    
     Z(7,1)=cte_s*sum(norms.A(7,1:3),2) ... 
        +cte_u*sum(abs(Am(3+3*m+1,4)) + abs(Am(3+3*m+1,3+m+1)) + abs(Am(3+3*m+1,3+2*m+1))+ abs(Am(3+3*m+1,3+3*m+1)))...
        +cte_u*sum(abs(Am(3+3*m+1:3+4*m,3+3*m+1)))...
        +params.L*abs(Am(3+3*m+1, 3+3*m+1:3+4*m))*((2*params.nu*a1vI(1:m)'+...
        2*params.nu*a1vI(3:m+2)' + 3*a1vIvI(1:m)' +3*a1vIvI(3:m+2)').*nu_power')...
        + z_inf(7);
    

    
   Z(1:7,2)=norms.A(1:7,1)*LAMBDA_S_1+norms.A(1:7,2)*LAMBDA_S_2+norms.A(1:7,3)*LAMBDA_S_4+...
        +norms.A(1:7,4)*LAMBDA_U_1+norms.A(1:7,5)*LAMBDA_U_2+norms.A(1:7,6)*LAMBDA_U_3+...
        +norms.A(1:7,7)*(LAMBDA_U_4+params.L*(nu+1/nu)*(8*params.nu+256*norms.a1));

   Z(1:7,3)=norms.A(1:7,7)*48*params.L*(nu+1/nu);

   disp('The last part of the bound Z is ')
    disp(mid(Z))

end
