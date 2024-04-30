function coeffs = getStBundleCoefficients(params)
    order = params.mfld.order; 
 %   eq = zeros(1, 4); 
 %   Df0 = JacSH_toMerge(0, params); 

    vectors.s = params.eigenvectors.s; 
    values.s = params.eigenvalues.s;
    vectors.u = params.eigenvectors.u; 
    values.u = params.eigenvalues.u;

   % v0 = vectors.s(:, 1)*params.scale; 
    lam1 = values.s(1); 
    lam2 = values.s(2); 

    Omega = diag([values.s(1), values.s(2), values.u(1), values.u(2)]);

    
    mflds.coeff.s = calc_proj_coeff(values.s, vectors.s, params); 
    A = DFQbundle(params, mflds); 
    Q0 = [vectors.s, vectors.u]; 


    %coeffs = zeros(4, order + 1, order + 1); 
    %coeffs(:, 1, 1) = v0; 
    
    % original coordinates 
    coeffs = zeros(4,2,order + 1, order + 1);
    coeffs(:,1,1,1) = vectors.s(:,1)*params.scale; 
    coeffs(:,2,1,1) = vectors.s(:,2)*params.scale;
    
    % eigenbasis coordinates 
    eigenbasis_coeffs = zeros(4,2,order + 1, order + 1); 
    eigenbasis_coeffs(:,1,1,1) = Q0^(-1)*vectors.s(:,1)*params.scale; 
    eigenbasis_coeffs(:,2,1,1) = Q0^(-1)*vectors.s(:,2)*params.scale;

    for alpha = 1:order
        for j = 0:alpha 
            
            i = alpha - j;
            %s_ij = starhatMat(A,coeffs, i,j); 
            disp(i)
            disp(j)
            
            % equation 20
            val = starhatMat(A,coeffs, i,j);
            % 4 x 2 vector
            s_ij = Q0^(-1)*val;
            %mat = ((i*lam1 + j*lam2 + growth_rate)*eye(4) ...
              %  - Df0)^(-1);
            % disp('s_ij')
            % disp(s_ij)

            % trying to implement equation 24
            mat1 = ((i*lam1 + j*lam2 + lam1)*eye(4) - Omega)^(-1);
            mat2 = ((i*lam1 + j*lam2 + lam2)*eye(4) - Omega)^(-1);

            coeff1 = mat1*s_ij(:,1);
            coeff2 = mat2*s_ij(:,2);

            eigenbasis_coeffs(:, :, i + 1, j + 1) = [coeff1, coeff2];
            coeffs(:,:,i+1, j+1) = [Q0*coeff1, Q0*coeff2]; 

          %  keyboard;
        end
    end








end
%     order = params.mfld.order;
% 
%     %order = 3;
%     equilibrium = zeros(1,4);
%     [vectors, values]= getJacEigs_toMerge(equilibrium, params);
% 
%     v0 = [vectors.u, vectors.s];
% 
%     v0 = v0*params.scale;
% 
%     %mflds.coeff.u=calc_proj_coeff(values.u,vectors.u,params);
%     mflds.coeff.s=calc_proj_coeff(values.s,vectors.s,params);
% 
%     %v0 = (0.0000001).*v0;
%     Omega = diag([values.u, values.s]);
% 
%     lam1 = values.s(1);
%     lam2 = values.s(2); 
% 
%     A = DFQbundle(params, mflds); 
% 
%     coeffs = zeros(4,4, order+1,order+1);
%     coeffs(:,:,1,1) = v0;
% 
%     % break up the computation to avoid recalculating the (0,0)
%     % coefficient. 
% 
%     i = 0;
%     for j = 1:order-i
%             s_ij = starhatMat(A,coeffs, i,j); 
%             tildes_ij = v0^(-1)*s_ij;
% 
%             % solve for the coefficients column by column
%             for k = 1:4 
%                 tildesk_ij = tildes_ij(:,k); 
%                 wk_ij = ((i*lam1+j*lam2+Omega(k,k))*eye(4) - Omega)^(-1)*tildesk_ij;
%                 coeffs(:,k,i+1,j+1) = v0*wk_ij;                 
%             end
% 
%     end
% 
%     for i = 1:order 
%         for j = 0:order-i
%             s_ij = starhatMat(A,coeffs, i,j); 
%             tildes_ij = v0^(-1)*s_ij;
% 
%             % solve for the coefficients column by column
%             for k = 1:4 
%                 tildesk_ij = tildes_ij(:,k); 
%                 wk_ij = ((i*lam1+j*lam2+Omega(k,k))*eye(4) - Omega)^(-1)*tildesk_ij;
%                 coeffs(:,k,i+1,j+1) = v0*wk_ij;                 
%             end
% 
%         end
%     end 
% end
