x = zeros(1, 4); 
params.lambda = 0; 
params.mu = 0.05; 
params.nu = 1.6; 
params.scale = 1/10;

order = 15; 
params.order = order; 
params.mfld.order = order; 

[eigenvectors, eigenvalues]= getJacEigs_toMerge(x, params);

params.eigenvalues.s = eigenvalues.s; 
params.eigenvectors.s = eigenvectors.s; 

st_mfld_coeffs = calc_proj_coeff(eigenvalues.s, eigenvectors.s, params);

lam_dot_alpha = lam_dot_alpha_mat(eigenvalues.s, order);
lam_dot_alpha = lam_dot_alpha(3: end, 3: end); 

d_dx_mfld_coeff = zeros(order + 1, order + 1, 4); 
for i = 1:4
    d_dx_mfld_coeff(:,:,i) = reshape(st_mfld_coeffs(:,:,i), [order + 1, order + 1]).*lam_dot_alpha; 
end

x_range = 0:.1:15; 
num = max(size(x_range)); 

d_dx_P = zeros(4, num); 
P = zeros(4, num); 
fcircP = zeros(4, num);
for i = 1:num 
    for j = 1:4
        d_dx_mat = exp(lam_dot_alpha.*x_range(i)).*reshape(d_dx_mfld_coeff(:,:,j), [order + 1, order + 1]);
        d_dx_P(j,i) = sum(sum(d_dx_mat)); 
        P_mat = exp(lam_dot_alpha.*x_range(i)).*st_mfld_coeffs(:,:,j); 
        P(j,i) = sum(sum(P_mat)); 
    end
    this_P = P(:,i);
    nonlin = params.nu*this_P(1)^2 - this_P(1)^3;
    this_f_circ_p = [this_P(2), ...
        this_P(3), ...
        this_P(4), ...
        -2*this_P(3) - (1 + params.mu)*this_P(1) + nonlin];
    fcircP(:, i) = this_f_circ_p; 
end 

diff = P - fcircP;

diff2 = d_dx_P - fcircP; 

delta = .1; 
gamma = P(1,:); 

num_diff = (gamma(2:end) - gamma(1:end-1))/delta;


figure; 
tiledlayout(4,1); 
nexttile; 
plot(x_range, real(d_dx_P(1,:))); 
title("d_dx_P")
nexttile; 
plot(x_range(1:end-1), real(num_diff))
title("num derivative")
nexttile
plot(x_range, real(fcircP(1,:)))
title("f circ P")
nexttile 
plot(x_range, real(P(1, :)))
title("P")

figure; 
tiledlayout(3,1); 
nexttile; 
plot(x_range, imag(d_dx_P(1,:))); 
nexttile; 
plot(x_range, imag(P(1,:))); 
nexttile;
plot(x_range, imag(fcircP(1,:)))
title('imaginary')
