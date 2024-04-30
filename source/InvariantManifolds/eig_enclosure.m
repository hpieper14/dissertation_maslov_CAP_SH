%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function constructs a radii poly for the eigenvalues and
% eigenvectors for the linearization of the SH equation at an equilibria
function poly=eig_enclosure(point,params,eval,evect,rstar)
    if evect(1)==0
        error('The first component is zero, try a different component')
    end
    
    ipoint=intval(point);
    iparams=intval(params);
    ieval=intval(eval);
    ievect=intval(evect);
    irstar=intval(rstar);
    
    DF=JacSH(ipoint(1),iparams(1),iparams(2));
    M=DF-ieval*intval(eye(4));
    
    bigMat=[-ievect,M(:,2:4)];
    A=bigMat^(-1);
    tildeZ0=0;
    Z2=2*max(abs(A),[],'all');
    tildeY0=max(abs(A*M*ievect), [],'all');
    
    % Now we construct eta -- we are just setting eta to be the largest
    % second derivative on the ball of radius rstar about the equilbria
    bigHess=intval(BigHessSH(ipoint(1)-irstar,iparams(2)));
    eta=[];
    % determines the component of the SH system for which we consider
    % the Hessian 
    for k=1:4
        % Now sum over the corresponding matrix 
        for i=1:4
            for j=1:4
                eta = [eta; abs(bigHess(i+4*(k-1), j))];
            end
        end
    end
    eta=max(eta);
    
     
    % Now we can construct the rest of the coefficients for the radii
    % polynomial
    
    % We need to make sure ||D^2f(c)|| 
    Y0=sup(tildeY0+eta*max(abs(A),[],'all')*max(abs(ievect))*irstar);
    Z0=sup(tildeZ0+eta*max(abs(A),[],'all')*irstar);
    
    % Now we search for a value at which the polynomial is negative
    poly=@(r)sup(Z2)*r^2-(1-Z0)*r+sup(Y0);
    
      
    p=[sup(Z2) -(1-Z0) sup(Y0)];
    val=roots(p);
    
    poly=min(val);
%     negative=false;
%     i=0;
%     while negative==false
%         val=intval(i*rad/1000);
%         if sup(poly(val))<0
%             val=sup(val);
%             negative=true;
%         end
%         i=i+1;
%         if i==1000
%             break 
%         end
%     end
%     
%     disp('Eigenvalue:')
%     disp(eval)
%     disp('Eigenvector:')
%     disp(evect)
%     
%     
%     if negative==true
%         fprintf('We found a negative value at r0=%d. \n',val)
%     else
%         disp('No negative value was found.')
%     end
%     poly=val;
end