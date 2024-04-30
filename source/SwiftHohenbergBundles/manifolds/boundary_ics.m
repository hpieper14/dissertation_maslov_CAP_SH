function ics = boundary_ics(coeff,p, theta, params) 
    theta=linspace(theta(1),theta(2), p);
    plotpoints=zeros(p,5);
    
     for k=1:p
        ps1s2 = zeros(4,1);
        s1=cos(theta(k));
        s2=sin(theta(k));
        for n=0:params.mfld.order
            for m=0:n
                point=reshape(coeff(n-m+1,m+1,:),[4,1]);
                ps1s2=ps1s2+point*(s1+1i*s2)^(n-m)*(s1-1i*s2)^m; 
            end
        end 
        plotpoints(k,:)=[real(ps1s2);theta(k)];
     end
     ics=plotpoints; 
     %plot3(plotpoints(:,1),plotpoints(:,2), plotpoints(:,4), 'o')
end