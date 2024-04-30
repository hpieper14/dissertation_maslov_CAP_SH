function sol=burke_knobloch_nf(params,phi,L)
%q_2=3/4-19*params.nu^2/18;
% I think this term is for the cubic quintic SH 

% this is taken from "Localized states in the generalized Swift-Hohenberg
% equation" by Burke and Knobloch

gamma3 = 38*params.nu^2/9 - 3;
x=-L:.01:L;
N=max(size(x));
y = 2*sqrt(2.*params.mu/gamma3).*sech(x.*sqrt(params.mu)./2).*cos(x+phi);
sol=zeros(N,2);
sol(:,1)=x.';
sol(:,2)=y.';
end
