function test = testBKNormalForm4dim()
    fourier.order = 10;
    params.nu = 1.2;
    params.mu = .1;
    vfParams = params;
    normalForm.branch = 0; 
    time = 10; 
    fourier.M = 1000;

    S = PulseSolution(fourier, vfParams, normalForm, time);
    
    % get data for normal form solution
    S = S.BKNormalForm4d_halfline(); 

    time = S.normalForm.time;
    sol = S.normalForm.sol;

    figure 
    tiledlayout(4,1)
    nexttile
    plot(time, sol(:, 1))
    nexttile 
    plot(time, sol(:, 2))
    nexttile 
    plot(time, sol(:, 3))
    nexttile 
    plot(time, sol(:, 4))
    title("Normal Form Solution")
    
    % trim time domain so normal form solution satisfies Neumann BC 
    S = trimNFSol_halfline(S); 

    time = S.normalForm.time;
    sol = S.normalForm.sol;

    figure 
    tiledlayout(4,1)
    nexttile
    plot(time, sol(:, 1))
    nexttile 
    plot(time, sol(:, 2))
    nexttile 
    plot(time, sol(:, 3))
    nexttile 
    plot(time, sol(:, 4))
    title("Normal Form Solution with Neumann BC")

    test = 0; 
end