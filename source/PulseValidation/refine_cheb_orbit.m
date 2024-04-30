function x = refine_cheb_orbit(x, mflds, params)
    disp('Refining Chebyshev coefficients via Newtons method.')
    k=0;
    m=params.cheb.order;
    while k < 150
        disp(['Newton iteration: ', num2str(k)])
                    
        fcn = FHomoclinic(x, mflds, params);
        DF = DF_homoclinic(x, params, mflds);
       
        fcn_vec=[fcn{1}, fcn{2}, fcn{3}, fcn{4}, fcn{5}, fcn{6}, fcn{7}];

        if vecnorm(fcn_vec)< params.tol 
            disp('Tolerance met.')
            break
        end
        
        if mod(k,5) == 0
            G2=FHomoclinic(x, mflds, params);
            disp('----------------------------------------------')
            disp(['Norm of F at iteration ',num2str(k),': ', num2str(vecnorm(fcn_vec))])
            disp(['Norm of DF^(-1)F at iteration ',num2str(k),': ', num2str(vecnorm(real(DF^(-1)*fcn_vec')))])
            disp(' ')
            disp('The value of F at the current x is: ')
            disp(['F1: ', num2str(G2{1})])
            disp(['F2: ', num2str(G2{2})])
            disp(['F3: ', num2str(G2{3})])
            
            disp('The max elements of (F4)_n,...,(F7)_n are: ')
            disp(['max(F4)_n: ', num2str(max(G2{4}))])
            disp(['max(F5)_n: ', num2str(max(G2{5}))])
            disp(['max(F6)_n: ', num2str(max(G2{6}))])
            disp(['max(F7)_n: ', num2str(max(G2{7}))])
            disp(' ')
            nozeros=abs(nonzeros(x.a4));
            max_ord=my_order(max(nozeros));
            min_ord=my_order(min(nozeros));
            disp(['Decay in Chebyshev coefficients F7: ', num2str(max_ord-min_ord)])
            disp('----------------------------------------------')
        end
        
        x_vec=[x.phi1, x.phi2, x.psi, x.a1, x.a2, x.a3, x.a4];
        realpart=real(DF^(-1)*fcn_vec');
        x_vec=x_vec-realpart';
        
        x.phi1=x_vec(1);
        x.phi2=x_vec(2); 
        x.psi=x_vec(3);
        x.a1=x_vec(3+1:3+m);
        x.a2=x_vec(3+m+1:3+2*m);
        x.a3=x_vec(3+2*m+1:3+3*m);
        x.a4=x_vec(3+3*m+1:3+4*m);
        
        k=k+1;
    end
    disp('Error after performing Newtons method: ')
    disp(vecnorm(real(DF^(-1)*fcn_vec')))
end

