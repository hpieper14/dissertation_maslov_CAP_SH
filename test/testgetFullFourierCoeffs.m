function test = testgetFullFourierCoeffs()
    order = 100;
    params.nu = 1.2;
    params.mu = .1;
    vfParams = params;
    nfBranch = 0; 
    time = 25; 

    S = PulseSolution(order, vfParams, nfBranch, time);
    S = BKNormalForm4dim(S);
    S = trimNFSol(S);
    S = getFullFourierCoeffs(S);
    
    fun = getFunctionFromFourierCoeffs(S);
 
    time = S.nfData.time; 
    sol = S.nfData.sol(:,1);
    figure
    plot(time, sol);
    hold on
    plot(fun(:,1),fun(:,2))
    legend('Sol via NF equation','Fourier approximation')

end