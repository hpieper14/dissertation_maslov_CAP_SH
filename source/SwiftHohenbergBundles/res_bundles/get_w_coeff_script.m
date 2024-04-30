x = zeros(1, 4); 
params.lambda = 0; 
params.mu = 0.05; 
params.nu = 1.6; 
params.scale = 1/3; 

order = 8; 

params.mfld.order = order; 

st_coeffs = getStBundleCoefficients(params);
v_coeff_pm = reshape(st_coeffs(1,:,:), [order + 1, order + 1]); 


[eigenvectors, eigenvalues] = getJacEigs_toMerge(0, params); 
st_mfld_coeffs = calc_proj_coeff(eigenvalues.s, eigenvectors.s, params);

mfld_deriv_coeff = st_mfld_deriv(params); 
v_coeff_dp = reshape(mfld_deriv_coeff(1,:,:), [order, order]); 
bottom = zeros(1, order); 
right = zeros(order + 1, 1); 
top = zeros(2, order + 1); 
left = zeros(order + 3, 2);
v_coeff_dp = [v_coeff_dp; bottom]; 
v_coeff_dp = [v_coeff_dp, right]; 
v_coeff_dp = [top; v_coeff_dp];
v_coeff_dp = [left, v_coeff_dp]; 



% pad v so it starts at index -2 
top = zeros(2, order + 1);
left_side = zeros(order + 3, 2);
v_coeff_pm = [top; v_coeff_pm]; 
v_coeff_pm = [left_side, v_coeff_pm];

w_coeff_pm = resonant_dprod_coeff(eigenvalues, v_coeff_pm, order - 3, - eigenvalues.s(1));
w_coeff_dp = resonant_dprod_coeff(eigenvalues, v_coeff_dp, order - 3, - eigenvalues.s(1));


lam_alpha_coeff = lam_dot_alpha_mat(eigenvalues.stable, order);
lam_alpha_coeff(3,3) = 1; 

prod_coeff_pm = w_coeff_pm.*lam_alpha_coeff; 
prod_coeff_dp = w_coeff_dp.*lam_alpha_coeff; 





% plot norms for all of these coefficients 

% plotpoints_w_pm=[];
% plotpoints_w_dp=[];
% plotpoints_v_pm=[];
% plotpoints_v_dp=[];
% 
% order = 6; 
% suborder=-2;
% while suborder<order
%     for i=-2:suborder
%         for j = suborder-i
%             if (j < order) && (i < order)
% 
% 
%             w_pm_norm = norm(w_coeff_pm(i+3, j+3)); 
%             w_dp_norm = norm(w_coeff_dp(i+3, j+3)); 
%             v_pm_norm = norm(v_coeff_pm(i+3, j+3)); 
%             v_dp_norm = norm(v_coeff_dp(i+3, j+3)); 
% 
% 
%             plotpoints_w_pm = [plotpoints_w_pm, [suborder; w_pm_norm]];
%             plotpoints_v_pm = [plotpoints_v_pm, [suborder; v_pm_norm]];
%             plotpoints_w_dp = [plotpoints_w_dp, [suborder; w_dp_norm]];
%             plotpoints_v_dp = [plotpoints_v_dp, [suborder; v_dp_norm]];
%             end
% 
%         end
%     end
%     suborder=suborder+1;
% end
% 
% figure 
% tiledlayout(2,1)
% nexttile 
% plot(plotpoints_v_pm(1,:), log(plotpoints_v_pm(2,:)),'o');
% xlabel('Coefficient order, $N = m + n$', Interpreter = 'latex')
% ylabel('$\log\left(\|p_{mn}\|\right)$', Interpreter = 'latex')
% title('Logarithmic Norms of Stable Manifold Coefficients $p_{mn}$', Interpreter = 'latex')


