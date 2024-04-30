% INTEGRATEDE Solves the normalized or non-normalized variational equation
% with initial condition ic on the interval [start, finish] and returns the
% pulse solution and the solution to the variational equation
%
%   [phisol, vsol] = integrateDE(C, ic,start,finish)  
%   [phisol, vsol] = C.integrateDE(ic,start,finish)  
function [phisol, vsol] =integrateDE(C, ic,start,finish)  
    % options for ODE45
    options=odeset('MaxStep',0.001, 'RelTol', 1e-13);

    if C.Euminus.normalize == 1
        [t,P]=ode45(@(t,P)C.nonautonODENormalized(t,P),[start finish],ic, options);
    else
        [t,P]=ode45(@(t,P)C.nonautonODE(t,P),[start finish],ic, options);
    end
    
    phisol=[t,P(:,1:4)];
    vsol=[t,P(:,5:8)];
end
