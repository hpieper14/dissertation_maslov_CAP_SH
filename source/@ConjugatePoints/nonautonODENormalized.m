% NONAUTONODENORMALIZED  Computes the normalized solution to the linearized first order system. 
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
% See also NONAUTONODE
function f1 = nonautonODENormalized(C, T, Y)
    params = C.vfParams; 
    lam=params.lambda;
    mu=params.mu;
    nu=params.nu;

    f1=Y(2);
    f2=Y(3);
    f3=Y(4);
    f4= -2.*Y(3) - (mu+1).*Y(1) + nu.*Y(1).^2 - Y(1).^3;
    
    B_x = [0,0,0,1;0,0,1,-2;lam-1-mu+2*nu.*Y(1)-3.*Y(1).^2,0,0,0;0,1,0,0];

    
    barw = [Y(5); Y(6); Y(7); Y(8)];
    
    vec = B_x*barw - barw.*dot(B_x*barw, barw);


    f1=[f1;f2;f3;f4;vec];
end
