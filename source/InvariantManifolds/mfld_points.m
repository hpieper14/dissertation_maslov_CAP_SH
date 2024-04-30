function plotpoints=mfld_points(coeff,params)
    p=30;
    
    % We plot them in polar coordinates 
    r=linspace(0,1,p);
    theta=linspace(0,2*pi, p);

    plotpoints=zeros(p,p,6);
    
    for j=1:p
        for k=1:p
            ps1s2 = zeros(4,1);
          
            s1=r(j)*cos(theta(k));
            s2=r(j)*sin(theta(k));
        for n=0:params.mfld.order
            for m=0:n
                point=reshape(coeff(n-m+1,m+1,:),[4,1]);
                ps1s2=ps1s2+point.*(s1+1i*s2)^(n-m).*(s1-1i*s2)^m;
                
            end
           if norm(ps1s2-real(ps1s2)) > 1e-10
               disp('There might be a problem. The imaginary part of the point is not negligible.')
           end
        end 
        plotpoints(j,k,:)=[real(ps1s2);r(j);theta(k)];
       
        end
    end
end