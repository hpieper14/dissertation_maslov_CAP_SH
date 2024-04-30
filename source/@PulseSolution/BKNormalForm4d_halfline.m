% BKNORMALFORM4D_HALFLINE  Computes the solution to the first order system
% from the normal form approximation to the pulse
% 
%   S = BKNormalForm4d_halfline(S)
%   S = S.BKNormalForm4d_halfline()
%
function S = BKNormalForm4d_halfline(S)
    r = S.vfParams.mu; 
    nu = S.vfParams.nu;
    gam = 38*nu^2/9-3;

    phi = S.normalForm.branch;
        
    M = S.fourier.M;
    dt = 2*pi/M; 
    t = dt*(1:M)';
    full_x=S.time*(t-pi)/pi;
    x = full_x(M/2:end);
    
    N=max(size(x));
    sol=zeros(N,4);

 
    sol(:,1) = 2*sqrt(2.*r./gam).*sech(x.*sqrt(r)./2).*cos(x+phi);
    sol(:,2) = - (2.*sin(phi + x).*(2.*r/gam)^(1/2))./cosh((r^(1/2).*x)./2)...
        - (r^(1/2).*sinh(r^(1/2).*x./2).*cos(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2).*x./2).^2;
    sol(:,3) = (r.*sinh(r^(1/2).*x./2).^2.*cos(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2).*x./2).^3 ...
        - (r.*cos(phi + x).*(2.*r./gam)^(1/2))./(2.*cosh(r^(1/2).*x./2)) ...
        - (2.*cos(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2)*x./2) ...
        + (2.*r^(1/2).*sinh(r^(1/2).*x./2).*sin(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2).*x./2).^2;
    sol(:,4) = (2.*sin(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2).*x./2) ...
        + (3.*r.*sin(phi + x).*(2.*r./gam)^(1/2))./(2.*cosh(r^(1/2).*x./2)) ...
        + (3.*r^(1/2).*sinh(r^(1/2).*x./2).*cos(phi + x).*(2.*r./gam)^(1/2))./cosh((r^(1/2).*x)./2).^2 ...
        + (5.*r^(3/2).*sinh(r^(1/2).*x./2).*cos(phi + x).*(2.*r./gam)^(1/2))./(4.*cosh(r^(1/2).*x./2).^2) ...
        - (3.*r.*sinh(r^(1/2).*x./2).^2.*sin(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2).*x./2).^3 ...
        - (3.*r^(3/2).*sinh(r^(1/2).*x./2).^3.*cos(phi + x).*(2.*r./gam)^(1/2))./(2.*cosh(r^(1/2).*x./2).^4);

    S.normalForm.time = x; 
    if S.vfParams.mu == .2
        S.normalForm.sol = 3*sol; 
    else 
        S.normalForm.sol = sol; 
    end
    
end
