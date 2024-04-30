function sol=BK_nf_4dim(params,phi,L, half)
    r = params.mu; 
    nu = params.nu;
    gam = 38*nu^2/9-3;
    
    if half == 1
        x=(-L:.01:0).';
    else
        x=(-L:.01:L).';
    end
    
        
    N=max(size(x));
    sol=zeros(N,5);
    
    sol(:,1)=x;
    sol(:,2) = 2*sqrt(2.*r./gam).*sech(x.*sqrt(r)./2).*cos(x+phi);
    sol(:,3) = - (2.*sin(phi + x).*(2.*r/gam)^(1/2))./cosh((r^(1/2).*x)./2)...
        - (r^(1/2).*sinh(r^(1/2).*x./2).*cos(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2).*x./2).^2;
    sol(:,4) = (r.*sinh(r^(1/2).*x./2).^2.*cos(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2).*x./2).^3 ...
        - (r.*cos(phi + x).*(2.*r./gam)^(1/2))./(2.*cosh(r^(1/2).*x./2)) ...
        - (2.*cos(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2)*x./2) ...
        + (2.*r^(1/2).*sinh(r^(1/2).*x./2).*sin(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2).*x./2).^2;
    sol(:,5) = (2.*sin(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2).*x./2) ...
        + (3.*r.*sin(phi + x).*(2.*r./gam)^(1/2))./(2.*cosh(r^(1/2).*x./2)) ...
        + (3.*r^(1/2).*sinh(r^(1/2).*x./2).*cos(phi + x).*(2.*r./gam)^(1/2))./cosh((r^(1/2).*x)./2).^2 ...
        + (5.*r^(3/2).*sinh(r^(1/2).*x./2).*cos(phi + x).*(2.*r./gam)^(1/2))./(4.*cosh(r^(1/2).*x./2).^2) ...
        - (3.*r.*sinh(r^(1/2).*x./2).^2.*sin(phi + x).*(2.*r./gam)^(1/2))./cosh(r^(1/2).*x./2).^3 ...
        - (3.*r^(3/2).*sinh(r^(1/2).*x./2).^3.*cos(phi + x).*(2.*r./gam)^(1/2))./(2.*cosh(r^(1/2).*x./2).^4);
    
end
