function satisfied = check_last_BC(rad, nu, x, params, mflds)
    m = params.cheb.order;
    N = params.mfld.order;
    
    irad = intval(rad);
    inu = intval(nu);    
    ia3 = intval(x.a3);
    
    idelta_s = intval(mflds.stable.error);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Interval around w3(1) %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    tail = ((1/inu)^m)/(intval(1)-intval(1)/inu);

    r_vec = irad./(inu.^(0:m-1)');
    ir_vec = infsup(-sup(r_vec),sup(r_vec));

    rad_w3 = (ia3(1)+ir_vec(1))+(intval(2)*sum(ia3(2:m)+ir_vec(2:m)'))+intval(2)*irad*infsup(-1,1)*tail;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interval around P3(theta) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ipsi = intval(x.psi);
    irho = intval(params.rho);
    ip3 = intval(mflds.stable.coeffs(:,:,3));
    
    rad_psi = infsup(inf(ipsi-irad),sup(ipsi+irad));
    itheta = [irho*cos(rad_psi);irho*sin(rad_psi)];
    
    % need to do this with a for loop 
    int_theta_powers = intval(zeros(N+1));
    
    for k = 0:N
        for j=0:N  
           int_theta_powers(k+1,j+1)=(itheta(1)+1i*itheta(2))^(k)*(itheta(1)-1i*itheta(2))^j;
        end
    end
    
    int_delta_s = infsup(-sup(idelta_s),sup(idelta_s));
    
    
    int_p3 = sum(sum(int_theta_powers.*ip3))+...
        +int_delta_s*((itheta(1)^(N+1)/(1-itheta(1)-1i*itheta(2)))*sum(((itheta(1)-1i*itheta(2))/(itheta(1)+1i*itheta(2))).^(0:N))+...
        +(1/(1-itheta(1)-1i*itheta(2)))*((itheta(1)-1i*itheta(2))^(N+1)/(1-itheta(1) + 1i*itheta(2))));

    % We verify that U2(1) and P2(theta) are both the same sign
    if (sup(rad_w3)<0 && sup(int_p3)<0) || (inf(rad_w3)>0 && inf(int_p3)>0)
        disp('Good to go! The boundary condition P3 = w3(1) is satisfied.')
        satisfied = 1;
    else
        disp('Stop! The last boundary condition is not satisfied.')
        satisfied = 0;
    end


end
