function sigma = getNewPulseMfldParams(oldMfldParams, params, desired_time, stable) 
    equilibrium = zeros(1,4);
    [~, values]= getJacEigs_toMerge(equilibrium, params);
    if stable == 1 
        Lambda = diag([values.s]);
    else 
        Lambda = diag([values.u]); 
    end
    sigma = exp(Lambda.*desired_time)*oldMfldParams;    
end
