function [x,params, mflds] = init_sol_val_tests()
m=450;


% these parameter values converge 

x.phi1=cos(5.755195195195195);
x.phi2=sin(5.755195195195195);
x.psi=3.669582882882883;
params.L = 3.922100000013819;
params.rho = 1 - .01;
params.scale = 3e-1;


x.a1=1:1:m;
x.a2=m+1:1:2*m;
x.a3=2*m+1:1:3*m;
x.a4=3*m+1:1:4*m;

params.mu=0.1;
params.nu=1.6;
params.lambda = 0;





params.cheb.order=m;
params.mfld.order=25;







params.tol=4e-16;

mflds=get_mflds(params);

end