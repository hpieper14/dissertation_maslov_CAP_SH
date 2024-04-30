function coeffs = computeFourierSpectrumNFSol(params, phi, type) 
% type corresponds to either using the normal form approximation or the
% integrated solution with the initial condition found via Chebyshev
    order = params.fourier.order;
    % tolerance for Newton's method
    tol = params.fourier.tol;
    L = params.fourier.L; 



    disp('----------------------------')
    disp('----------------------------')
    disp('Computing Fourier coeff for solution.')
 

    if type == 0
        sol=BK_nf_4dim(params,phi,L, 0);
    elseif (type == 1) && (phi == 0)
        normalize = 1;
        u_ic = [-0.190343484893489 ; -0.128234936941383 ; 0.193875175006471 ; 0.168311653800752];
        params.L = 3.922100000013819;
        [sol, frame] = generateEuFrame(0,params.fourier.L, params, u_ic);

        
    elseif (type == 1) && (phi == pi)
        normalize = 1;
        u_ic = [-0.190343484893489 ; -0.128234936941383 ; 0.193875175006471 ; 0.168311653800752];
        params.L = 3.922100000013819;
        [sol, frame] = generateEuFrame(0,params.fourier.L, params, u_ic);
        sol(:,2) = - sol(:,2);
    end
    




    % determine where to truncate the solution -- i.e. we want the first
    % derivative to be 0 and the function value to be 0

    index = find(diff(sign(sol(:,3)))); 
    pulse2 = sol(index:end-index,2);
    t2 = sol(index:end-index,1);
    L = -t2(1);


    test_coeffs = getFullFourierCoeffs(pulse2, t2, order, L);
    side = 2;


    disp('-----------------------------------------------')
    disp('Performing Newtons method to refine the coefficients.')




    betternfcoeffs=NewtonFourierNDBC(test_coeffs,order, params, tol);
    betternff=fourierODE(betternfcoeffs, order,params,side);



    disp('These coefficients should be a zero of F.')
    disp('Max elt F(a)=')
    disp(max(betternff))

    format long

    disp('The nonnegative eigenvalues of DF are:')
    DF=DFFourierNDBC(betternfcoeffs, order,params) ;
    D=eig(DF);
    [i,j] = find(D>0);
    if min(size(D(i,j))) == 1
        disp(D(i,j));
    else 
        A = D(i,j);
        disp(A(:,1))
    end
    

    format short

    T=-params.fourier.L:.05:params.fourier.L;
    f = 0;
    for n = -order:order
        f = f+real(betternfcoeffs(n+order+1)*exp(1i*n*pi.*T./params.fourier.L));
    end
  
    figure 
    plot(t2, pulse2)
    hold on
    plot(T,f)
    legend('Sol via NF','solution after Newton')
    title('Fourier approximation of solution after Newtons method')
    coeffs = betternfcoeffs;
end
