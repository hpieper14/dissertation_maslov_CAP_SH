%------------------------------------------------------------------------
% initialize parameters 
%------------------------------------------------------------------------

params.lambda = 0; 
params.mu = 0.1;
params.nu = 1.2; 

u_ic =  [-0.190343484893489 ; -0.128234936941383 ; 0.193875175006471 ; 0.168311653800752];
params.L = 3.922100000013819;

[vectors, values]= getBinfEigs(params);
v1=real(vectors.u(:,1));
v2=imag(vectors.u(:,1));
params.normalizeBasis = 0;

LB = 15;

%-----------------------------------------------------------------------
% Compute Cheb Coeffs 
%------------------------------------------------------------------------
withError = 0;
[full_phi1, frame] = generateEuFrame(withError, -10, params, u_ic);
sol_vec = full_phi1(:,2:5);
time_vec = full_phi1(:,1);

M = 100;

% must call chebfun on a column vector
x1=chebfun(sol_vec(:,1),'trunc', M);
x2=chebfun(sol_vec(:,2),'trunc', M);
x3=chebfun(sol_vec(:,3),'trunc', M);
x4=chebfun(sol_vec(:,4),'trunc', M);

a1=chebcoeffs(x1);
a2=chebcoeffs(x2);
a3=chebcoeffs(x3);
a4=chebcoeffs(x4);

cheb_coeff=[a1, a2, a3, a4];


figure
tiledlayout(4,1)
nexttile
plot(x1)
nexttile
plot(x2)
nexttile
plot(x3)
nexttile
plot(x4)
title('Cheb Solution (Pre-Newton Method)')


figure
tiledlayout(4,1)
nexttile
plot(time_vec, sol_vec(:,1))
nexttile
plot(time_vec, sol_vec(:,2))
nexttile
plot(time_vec, sol_vec(:,3))
nexttile
plot(time_vec, sol_vec(:,4))
title('Solution via Integrating')
    