function branch = pseudoArclengthCont(a,params, step)
    side = 2;
    order = params.fourier.order; 
    if size(a,1) == 1
        DFmu = -a';
        a = a';
    else 
        DFmu = -a;
    end
    
    A = [DFFourierNDBC(a,order, params), DFmu];
    X = null(A);
    if (size(X,2) == 1) == 0 
        X = X(:,2);
    end
    dotX=(X)./vecnorm(X);
 
    predictor = [a; params.mu] + step.*dotX;
    
    X = predictor;
    close_enough = false; 
    i = 0;
    while close_enough == false 
        %disp(['iteration:', num2str(i)])
        params.mu = X(end);
        a = X(1:end-1);
        E = dot(X - predictor, dotX); 
        F = vertcat(E, fourierODE(a,order,params, side)');
        if abs(vecnorm(F))<params.fourier.tol
            close_enough=true;
            disp('Tolerance met')
            break;
        end
        DFmu = -a;
        DFN = [DFFourierNDBC(a,order, params), DFmu];
        DF = [dotX'; DFN];
        
        X = X-((DF^(-1)*F));
        i = i+1;
    end
    branch = X; 
end
