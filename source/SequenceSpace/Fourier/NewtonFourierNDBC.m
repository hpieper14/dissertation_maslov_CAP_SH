function zero=NewtonFourierNDBC(s,order,params,tol)
    side = 2;
    k=0;
    while k < 350
        fcn=fourierODE(s, order,params,side);
        DF=DFFourierNDBC(s,order,params);
        if vecnorm(fcn)< tol 
            disp('tolerance met')
            break
        end
        s=s-(DF^(-1)*fcn')';
        k=k+1;
    end
    zero=s;
%    disp('Approx zero:')
%    disp(s)
    disp("used iterations: ")
    disp(k)
    disp('Error:')
    disp(vecnorm(fcn))   
end