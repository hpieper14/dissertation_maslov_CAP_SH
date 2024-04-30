 function [branch, norms] = makeSolBranch(step, Niter, orig_params, orig_sol)
    mu_init = orig_params.mu; 
    
    
    if size(orig_sol, 1) == 1
        orig_sol = orig_sol';
    end
    
    
    %N = max(size(mu_range));

    branch = [orig_sol; mu_init; vecnorm(orig_sol)];
    sol = getFunctionFromFourierCoeffs(orig_sol, orig_params.fourier.L, orig_params.fourier.order);
    fcn = sol(:,2);
    norms = [vecnorm(fcn)];
    a = orig_sol;
    params = orig_params;
    for i = 1:Niter
        newzero = pseudoArclengthCont(a,params, step);
        
        params.mu = newzero(end); 
        a = newzero(1:end-1); 
        sol = getFunctionFromFourierCoeffs(a, params.fourier.L, params.fourier.order);
        fcn = sol(:,2);
        norms = [norms, vecnorm(fcn)];
        newzero = [newzero; vecnorm(fcn)];
        
        branch = [branch, newzero];
        
        if mod(i, 10) == 0
            figure 
            plot(sol(:,1), sol(:,2))
            xlabel('time')
            ylabel('$\|u\|_2$', 'Interpreter', 'latex')
            title(['\mu =', num2str(branch(end-1, end))], 'Interpreter', 'latex');
        end
        
            
    end
    
    
    
    all_mu = branch(end-1, :); 
    format_data = sortrows([all_mu; norms]', 1);
    figure 
    plot(format_data(:,1),format_data(:,2),'o');
    xlabel('$\mu$','Interpreter', 'Latex');
    ylabel('$\|u\|$', 'Interpreter', 'Latex');
    
    
end

    