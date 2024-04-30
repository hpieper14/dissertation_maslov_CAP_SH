function pulses = TestPulseIntegrationMethods()
    
    % initialize parameters
    params.mu=.1;
    params.nu=1.6;
    params.lambda = 0;
    params.scale = 3e-1;
    params.L = 40;
    
    u_ic = [-0.000753881878776;   0.000675254086116;   0.001002574703915;  -0.000394228947926];

    options=odeset('MaxStep',0.001, 'RelTol', 1e-8);
    [t,P]=ode45(@(t,P)SH_comp(t,P,params),[-params.L params.L],u_ic, options);
    x = zeros(1,4);
    [vectors, values]= get_eigs(x, params);
    v1=real(vectors.u(:,1));

    ic = [u_ic; v1];
    normalize = 0;
    [full_phi1, basis_1] = integrateDE(ic,-params.L,params.L,params, normalize);
    
    figure 
    plot(t,P(:,1))
    title('Pulse solution via numerical integration of SH eq')
    
    figure
    plot(full_phi1(:,1), full_phi1(:,2))
    title('Pulse solution via numerical integration of coupled ODE (pulse + basis sol)')

end
