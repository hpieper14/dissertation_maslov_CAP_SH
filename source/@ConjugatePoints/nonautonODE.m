% NONAUTONODE  Computes the solution to the linearized first order system. 
% Since this is a nonautonomous system, we compute the pulse solution and 
% the solution to the variational equation simultaneously. 
%   f1 = nonautonODE(C, T, Y)
%   f1 = C.nonautonODE(T, Y)
% 
% returns a matrix of function values with rows [\varphi_1, \dots,
% \varphi_4, u_1, u_4] where \varphi_i is the ith component of the pulse
% solution and u_j is the jth component of the solution to the variational
% equation g'\circ\varphi
%
% See also NONAUTONODENORMALIZED
function f1 = nonautonODE(C, T, Y)
    params = C.vfParams; 
    lam=params.lambda;
    mu=params.mu;
    nu=params.nu;

    f1=Y(2);
    f2=Y(3);
    f3=Y(4);
    f4= -2.*Y(3) - (mu+1).*Y(1) + nu.*Y(1).^2 - Y(1).^3;

    g1 = Y(8);
    g2 = Y(7) - 2.*Y(8);
    g3=(lam - 1 - mu + 2*nu.*Y(1) - 3.*Y(1).^2).*Y(5);
    g4=Y(6);

    f1=[f1;f2;f3;f4;g1;g2;g3;g4];
end
