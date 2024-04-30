%[x,params, mflds] = init_sol_val_tests();

params.L = 3.922100000013819;
params.mu = .1;
params.nu = 1.6;
params.lambda = 0;

% u_ic = [-0.190343484893489 ; -0.128234936941383 ; 0.193875175006471 ; 0.168311653800752];
% 
% L = 50;
% 
% [V, D]=asym_un_vecs(params);
% v1=real(V(:,1));
% 
% normalize = 1;
% 
% [full_phi1, EuBasis] = generateEuFrame(L, params, u_ic, normalize);
% 
% time=full_phi1(:,1);
% full_phi1(:, 2:5) =full_phi1(:,2:5);
% 
% figure
% plot(time, full_phi1(:,2))
% title('Pulse')
% 
% 
% step=diff(time);
% u = full_phi1(:,2);
% uprime=full_phi1(:,3);
% u2prime=full_phi1(:,4);
% u3prime=full_phi1(:,5);
% u4prime=[diff(u3prime)./step; 0];
% 
% coupled_test_vec= -u4prime - 2.*u2prime - (1+params.mu).*u + params.nu.*u.^2 - u.^3;
% figure
% plot(time, coupled_test_vec)
% title('Plug coupled solution back into ODE -- should be 0')
% 
% half = 0;
% sol0=BK_nf_4dim(params,0,L, half);
% time = sol0(:,1);
% 
% 
% step=diff(time);
% u = sol0(:,2);
% uprime=sol0(:,3);
% u2prime=sol0(:,4);
% u3prime=sol0(:,5);
% u4prime=[diff(u3prime)./step; 0];
% 

params.fourier.order = 150;
params.fourier.L = 100; 
params.fourier.tol = 1e-12;
coeffs = computeFourierSpectrumNFSol(params,pi, 1);

sol = getFunctionFromFourierCoeffs(coeffs, params.fourier.L, params.fourier.order); 

step=diff(sol(:,1));
phi = sol(:,2);
phix = [diff(phi)./step; 0];
phixx = [diff(phix)./step; 0];
phixxx = [diff(phixx)./step; 0];
phixxxx = [diff(phixxx)./step; 0];

test_vec = -phixxxx - 2.*phixx - (1+params.mu).*phi + params.nu.*phi.^2 - phi.^3;
figure
plot(sol(:,1), test_vec)
title('Plug phi = pi back into ODE -- should be 0')


coeffs = computeFourierSpectrumNFSol(params,0, 1);

sol = getFunctionFromFourierCoeffs(coeffs, params.fourier.L, params.fourier.order); 

step=diff(sol(:,1));
phi = sol(:,2);
phix = [diff(phi)./step; 0];
phixx = [diff(phix)./step; 0];
phixxx = [diff(phixx)./step; 0];
phixxxx = [diff(phixxx)./step; 0];

test_vec = -phixxxx - 2.*phixx - (1+params.mu).*phi + params.nu.*phi.^2 - phi.^3;
figure
plot(sol(:,1), test_vec)
title('Plug phi = 0 back into ODE -- should be 0')

% [ut,ux]=ode45(@(st,sx) SH_comp(st,sx, params), [0,L1], u_ic, options);
% [st,sx]=ode45(@(st,sx) SH_comp(st,sx, params), [L1, 0], s_ic, options);
% 
% time=[-st;ut];
% u = [ux(:,1);flip(sx(:,1))];
% uprime=[ux(:,2); flip(sx(:,2))];
% u2prime=[ux(:,3); flip(sx(:,3))];
% u3prime=[ux(:,4); flip(sx(:,4))];
% u4prime=[0;diff([ux(:,4); flip(sx(:,4))])];
% 
% f = -mu.*u + nu.*u.^2 - u.^3;
% 
% test_vec1=-u4prime-2.*u2prime-u + f;

