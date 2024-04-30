% This function solves the ODE for the pulse and the basis solution for the
% unstable subspace

function [phisol, vsol] =integrateDE(ic,start,finish,params)  
    options=odeset('MaxStep',0.001, 'RelTol', 1e-10);
    if params.normalizeBasis == 0
        [t,P]=ode45(@(t,P)bigPhi(t,P,params),[start finish],ic, options);
    else
        [t,P]=ode45(@(t,P)bigPhi_normalized(t,P,params),[start finish],ic, options);
    end
    
    phisol=[t,P(:,1:4)];
    vsol=[t,P(:,5:8)];
end