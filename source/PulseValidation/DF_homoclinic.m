function DF = DF_homoclinic(x, params, mflds)
    m = params.cheb.order;
    DF = zeros(4*m + 3, 4*m + 3);   
    
    %%% Quantities in the derivatives of manifold parameterizations 
    
    phi=[x.phi1, x.phi2];
    theta = [params.rho*cos(x.psi), params.rho*sin(x.psi)];
    
    N = params.mfld.order;

    theta_pow1=zeros(N+1);
    theta_pow2=zeros(N+1);
    phi_pow1=zeros(N+1);
    phi_pow2=zeros(N+1);
    
    for k = 0:N
        for j=0:N
            phi1=k*(phi(1)+1i*phi(2))^(k-1)*(phi(1)-1i*phi(2))^j;
            phi2=j*(phi(1)+1i*phi(2))^k*(phi(1)-1i*phi(2))^(j-1);
            
            psi1=k*(theta(1)+1i*theta(2))^(k-1)*(theta(1)-1i*theta(2))^j;
            psi2=j*(theta(1)+1i*theta(2))^k*(theta(1)-1i*theta(2))^(j-1);
            
            theta_pow1(k+1,j+1)=psi1;
            theta_pow2(k+1,j+1)=psi2;
            
            phi_pow1(k+1,j+1)=phi1;
            phi_pow2(k+1,j+1)=phi2;
        end
    end
    
    % vectorized notation takes a complex conjugate and switches a sign? 
    %theta_pow1=((0:N).*((theta(1)+1i*theta(2)).^[0 (0:N-1)]))'*((theta(1)+1i*theta(2)).^(0:N));
    %theta_pow2=((theta(1)+1i*theta(2)).^(0:N))'*((0:N).*((theta(1)+1i*theta(2)).^[0 (0:N-1)]));

    %%% Convolutions in the derivative

    [a1a1]=chebstar2(x.a1, x.a1,2*m);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % UPPER LEFT 3x3 MATRIX %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    p1 = mflds.stable.coeffs(:,:,1);
    p2 = mflds.stable.coeffs(:,:,2);
    p4 = mflds.stable.coeffs(:,:,4);
    
    aprime = -params.rho*sin(x.psi)+1i*params.rho*cos(x.psi);
    bprime = -params.rho*sin(x.psi)-1i*params.rho*cos(x.psi);
    
    DF(1,3) = -aprime*sum(sum(theta_pow1.*p1))-bprime*sum(sum(theta_pow2.*p1));
    DF(2,3) = -aprime*sum(sum(theta_pow1.*p2))-bprime*sum(sum(theta_pow2.*p2));
    DF(3,3) = -aprime*sum(sum(theta_pow1.*p4))-bprime*sum(sum(theta_pow2.*p4));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOWER LEFT 4mx3 MATRIX %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    q1 = mflds.unstable.coeffs(:,:,1);
    q2 = mflds.unstable.coeffs(:,:,2);
    q3 = mflds.unstable.coeffs(:,:,3);
    q4 = mflds.unstable.coeffs(:,:,4);
    
    DF(4,1) = -sum(sum(phi_pow1.*q1)) - sum(sum(phi_pow2.*q1));
    DF(4,2) = -sum(sum(phi_pow1.*q1))*(1i)- sum(sum(phi_pow2.*q1))*(-1i);
                
    DF(3+m+1,1) = -sum(sum(phi_pow1.*q2))- sum(sum(phi_pow2.*q2));
    DF(3+m+1,2) = -sum(sum(phi_pow1.*q2))*(1i) - sum(sum(phi_pow2.*q2))*(-1i);
                
    DF(3+2*m+1,1) = -sum(sum(phi_pow1.*q3))- sum(sum(phi_pow2.*q3));
    DF(3+2*m+1,2) = -sum(sum(phi_pow1.*q3))*(1i)- sum(sum(phi_pow2.*q3))*(-1i);
                
    DF(3+3*m+1,1) = -sum(sum(phi_pow1.*q4)) - sum(sum(phi_pow2.*q4));
    DF(3+3*m+1,2) = -sum(sum(phi_pow1.*q4))*(1i) - sum(sum(phi_pow2.*q4))*(-1i);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UPPER RIGHT 3x4m MATRIX %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DF(1,3+1:3+m) = [1, 2*ones(1,m-1)];
    DF(2,3+m+1:3+2*m) = [1, 2*ones(1,m-1)];
    DF(3,3+3*m+1:3+4*m) = [1, 2*ones(1,m-1)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOWER RIGHT 4mx4m MATRIX %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pm_ones = (-1).^(1:m-1);
    
    %%%%
    %f4%
    %%%%
    
    % d(f4)_0 da1
    DF(3+1,3+1:3+m) = [1, 2*pm_ones]; 
    
    % d(f4)_k d(a1)_k, k \geq 1 (column vector)
    for i = 2:m
        DF(3+i, 3+i) = 2*(i-1);
    end
    
    % df4 da2
    for k = 1:m-1
        for l = 0:m-1
            if l==k-1
                DF(3+1+k, 3+m+l+1) = -params.L;
            elseif l==k+1
                DF(3+1+k, 3+m+l+1) = +params.L;
            end
        end
    end
    
    
    %%%%
    %f5%
    %%%%
    
    
    % d(f5)_0 da2
    DF(3+1+m,3+1+m:3+2*m) = [1, 2.*pm_ones]; 
    
    % d(f5)_k d(a2)_k, k \geq 1 (column vector)
    for i = 2:m 
        DF(3+m+i, 3+m+i) = 2*(i-1);
        
    end
    
    % df5 da3 
  for k = 1:m-1
        for l = 0:m-1
            if l==k-1
                DF(3+1+m+k, 3+2*m+l+1) = -params.L;
            elseif l==k+1
                DF(3+1+m+k, 3+2*m+l+1) = +params.L;
            end
        end
   end
  
    %%%%
    %f6%
    %%%%
    
   
    
    % d(f6)_0 da3
    DF(3+1+2*m,3+1+2*m:3+3*m) = [1, 2.*pm_ones]; 
    
    % d(f6)_k d(a3)_k, k \geq 1 (column vector)
    for i = 2:m
        DF(3+2*m+i, 3+2*m+i) = 2*(i-1);
    end
   
  % df6 da4
  for k = 1:m-1
        for l = 0:m-1
            if l==k-1
                DF(3+1+2*m+k, 3+3*m+l+1) = -params.L;
            elseif l==k+1
                DF(3+1+2*m+k, 3+3*m+l+1) = +params.L;
            end
        end
   end
   
    %%%%
    %f7%
    %%%%
    

    % df7_k da1_l, k >0 

    x.a1=[x.a1, zeros(1,size(x.a1,2))];

    for k = 1:m-1
        for l = 0:m-1
            if l == 0
                DF(4+3*m+k, 4+l) = -params.L*(2*params.nu*(x.a1(abs(k-1-l)+1)) ...
                                    - 3*(a1a1(abs(k-1-l)+1)) ...
                                    - 2*params.nu*(x.a1(abs(k+1-l)+1))...
                                    + 3*(a1a1(abs(k+1-l)+1)));
            else
                DF(4+3*m+k, 4+l) = -params.L*(2*params.nu*(x.a1(abs(k-1-l)+1)+x.a1(abs(k-1+l)+1)) ...
                                    -3*(a1a1(abs(k-1-l)+1)+a1a1(abs(k-1+l)+1)) ...
                                    - 2*params.nu*(x.a1(abs(k+1-l)+1)+x.a1(abs(k+1+l)+1))...
                                    + 3*(a1a1(abs(k+1-l)+1)+a1a1(abs(k+1+l)+1))); 
            end

            if l==k-1
                DF(4+3*m+k, 4+l) = DF(4+3*m+k, 4+l) - params.L*(-1-params.mu);
            elseif l==k+1
                DF(4+3*m+k, 4+l) = DF(4+3*m+k, 4+l) + params.L*(-1-params.mu);
            end
        end
    end


        % d(f7)_0 da4 
        DF(3+3*m+1, 3+3*m+1:3+4*m) = [1, 2.*pm_ones];

        %d(f7)_k d(a4)_k
        for i = 2:m 
            DF(3+3*m+i, 3+3*m+i) = 2*(i-1);
        end

       % d(f7)_k d(a3)_l
       for k = 1:m-1
            for l = 0:m-1
                if l==k-1
                    DF(3+1+3*m+k, 3+2*m+l+1) = 2*params.L;
                elseif l==k+1
                    DF(3+1+3*m+k, 3+2*m+l+1) = -2*params.L; 
                end
            end
       end

end
