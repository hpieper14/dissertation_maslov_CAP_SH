%clear all 

%-------------------------------------------------------------------------
% first we calculate the size of the eigenvector/eigenvalue enclosures
point=[0,0,0,0];
params=[0.05,1.6]; % in the order [mu,nu']
Df0=JacSH(0,params(1),params(2));
[V1,D1]=eigs(Df0);
[d, ind]=sort(real(diag(D1)));
D=D1(ind,ind);
V=V1(:,ind);

rstar=1e-15;

error=zeros(1,4);
for i=1:4
    error(i)=eig_enclosure(point,params,D(i,i),V(:,i),rstar);
end
error=max(error);

%----------------------------------------------------------------------
% Now we scale the eigenvectors in accordance with Theorem 10.5.1 so the
% last component is on the order of machine precision.
tau=1/7;
Vscale=V*tau;

% Separate the eigenvalues and associated vectors into those with positive
% and negative real part
uneigs=[D(1,1),D(2,2)];
unvec=Vscale(:,1:2);
stabeigs=[D(3,3),D(4,4)];
stabvec=Vscale(:,3:4);

% -----------------------------------------------------------------------
% Now we calculate the coefficients of the parameterization for the stable
% and unstable manifold up to a desired order. 
order=16;

% unstable
disp('Calculating the coefficients for the unstable manifold.')
uncoeff=calc_proj_coeff(order,uneigs,unvec,params);

% stable 
disp('Calculating the coefficients for the stable manifold.')
stabcoeff=calc_proj_coeff(order,stabeigs,stabvec,params);


k=6;
mat=zeros(k^2,4);
K = zeros(k^2, 2);
for i = 1:k 
    for j = 1:k
        row=(i-1)*k + j;
        mat(row,:)=stabcoeff(i,j,:);
        K(row,:)=[i-1,j-1];
    end  
end

figure
tiledlayout(2,2) 

nexttile
plot_coeff(uncoeff,order);
title('Coefficient norms for unstable parameterization')

nexttile
plot_coeff(stabcoeff,order);
title('Coefficient norms for stable parameterization')
 
nexttile
plot_manifold(uncoeff,order);
title('Unstable Manifold')

nexttile
plot_manifold(stabcoeff,order);
title('Stable Manifold')

%--------------------------------------------------------------------------
% Now we apply Lemma 10.4.1 to validate the parameterization we computed 


 unpoly=radiipoly(params,uncoeff,V,D,order);
 unstable_bound = min(unpoly(find(unpoly > 0)));

 stpoly=radiipoly(params,stabcoeff,V,D,order);
 stable_bound = min(stpoly(find(stpoly > 0)));
%  
  disp('The error on the parameterization for the unstable manifold is: ')
  disp(unstable_bound)
  disp('The error on the parameterization for the stable manifold is: ')
 disp(stable_bound)

%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function validates the parameterization of the invariant manifold 
% Inputs: params - 2x1 vector in the form [mu, nu']
%         coeff - multidimensional array of size 3*order+1 x 3*order+1 x 4.
%         This could correspond to the stable or unstable manifold
%         Q - matrix of unscaled eigenvectors 
%         Lambda - diagonal matrix of the eigenvalues 
%         order - order of the approximation - N=m+n
%         rad - the radius on which to search for a negative value
% Output: poly - value at which the sup of the polynomial evaluated with
%         interval arithmetric is negative
function val=radiipoly(params,coeff, Q,Lambda,order)
    % Calculate the coefficient K 
    maxKmn=((order+1)*abs(real(Lambda(1))) - abs(Lambda(1)))^(-1);
 
    K_N = max(abs(Q),[],'all')*max(abs(Q^(-1)),[],'all')*maxKmn;
    
    % extract the coefficients a_{mn}
    
    a=coeff(1:order+1,1:order+1,1);
    
    A=zeros(order+1, 2*order);
    B=zeros(2*order, order+1);
    C=zeros(2*order, 2*order);
    
    bara=[a, A ; B, C];
    %bara=coeff(1:3*order+1, 1:3*order+1,1);
    
    disp('Calculating Y0.')

    %%%%%%
    % Y0 %
    %%%%%%

    quadsum=0;
    suborder=order+1;
    while suborder < 2*order+1
        for i=0:suborder
            j=suborder-i;
            astarij=starhat(bara, bara, i,j);
            quadsum=quadsum+abs(astarij);
        end
        suborder=suborder+1;
    end
   
    cubsum=0;
    suborder=order+1;
    while suborder<3*order+1
        for i=0:suborder
            j=suborder-i;
            tripstara=tripstarhat(bara,bara,bara,i,j);
            cubsum=cubsum+abs(tripstara); 
        end
        suborder=suborder+1;
    end
    
    Y0 = K_N*(params(2)*quadsum+cubsum);
    disp(Y0)
    disp('Calculating Z1.')
    
    %%%%%%
    % Z1 %
    %%%%%%
    
    lowerquadsum=0;
    suborder=1;
    while suborder < order+1
        for i=0:suborder
            j=suborder-i;
            astarij=starhat(bara, bara, i,j);
            lowerquadsum=lowerquadsum+abs(astarij);
        end
        suborder=suborder+1;
    end
    
    linsum=0;
    suborder=1;
    while suborder<order+1
        for i=0:suborder
            j=suborder-i;
            linsum=linsum+abs(bara(i+1,j+1));
        end
        suborder=suborder+1;
    end
    
    Z1= K_N*(2*params(2)*linsum + 3*lowerquadsum);
    disp(Z1);
    disp('Calculating Z2.')
    %%%%%%
    % Z2 %
    %%%%%%
    
    Z2=@(r)K_N*(6*linsum+2*params(2)+3*r);
    disp(Z2)
    
    b=K_N*3;
    a=K_N*(6*linsum + 2*params(2));
  
    disp(a)
    disp(b)
    % this section of code builds an interval on which the sup of the polynomial is negative 
    poly=@(r)(Z1+Z2(r)*r)*r+Y0-r;
    
    %figure(3)
    %fplot(@(r) poly(r));
    
    p=[b, a, Z1-1, Y0];
    
    val=roots(p);
    
    % 
    % 
    % negative=false;
    % i=0;
    % while negative==false
    % 
    % 
    %     val=intval(i*rad/1000);
    %     if sup(poly(val))<0
    %         val=sup(val);
    %         negative=true;
    %     end
    %     i=i+1;
    %     if i==1000
    %         break 
    %     end
    % end
    % if negative==true
    %     fprintf('We found a negative value at r0=%d. \n',val)
    % else
    %     disp('No negative value was found.')
    % end
    poly=val;
end


% This function evaluates the function P on a grid of points and plots them
% Inputs: coeff - corresponds to the coefficients computed for either the
%                 stable or unstable manifold
%         order - order of the parameterization 
% Outputs: plots = 0
%          also generates a 3D plot. In this case, we omit the third
%          component (seemed to generate the best plot)
function plots=plot_manifold(coeff,order)
    p=10;
    theta=linspace(0,2*pi,p);
    r=linspace(0,1,p);
    plotpoints=zeros(p,p,4);
    
    for j=1:p
        for k=1:p
        ps1s2 = zeros(4,1);
        s1=r(j)*cos(theta(k));
        s2=r(j)*sin(theta(k));
        for n=0:order
            for m=0:n
                %disp(n)
                %disp(m)
                point=reshape(coeff(n-m+1,m+1,:),[4,1]);
                ps1s2=ps1s2+point*(s1+1i*s2)^(n-m)*(s1-1i*s2)^(m);
                
            end
        end 
        plotpoints(j,k,:)=real(ps1s2);
        end
    end
    surf(plotpoints(:,:,1),plotpoints(:,:,2),plotpoints(:,:,4), plotpoints(:,:,3));
    xlabel('x1')
    ylabel('x2')
    zlabel('x4')
    plots=0; 
end 

% this function plots the log of the norms of the coefficients against
% their order. This is done to check for decay. 
% Inputs: coeff - corresponds to the coefficients computed for either the
%                 stable or unstable manifold
%         order - order of the parameterization 
% Outputs: plots = 0
%          also generates a 2D plot
function plots=plot_coeff(coeff,order)
    plotpoints=[];
    suborder=1;
    while suborder<order+1
        for i=0:suborder
            for j = suborder-i
                vec=zeros(1,4);
                for k=1:4
                    vec(k)=coeff(i+1,j+1,k);
                end
                normpoint=norm(vec);
                point=[suborder;normpoint];
                plotpoints=[plotpoints,point];
            end
        end
        suborder=suborder+1;
    end
    
  plot(plotpoints(1,:), log(plotpoints(2,:)),'o');
  xlabel('Order of the coefficients')
  ylabel('Log norm of the coefficients')
  plots=0;
end

% This function recursively calculates the coefficients of the
% parameterization for the invariant manifold 
% Inputs: order - order of the desired parameterization 
%         eigenvalues - either the stable or unstable pair of eigenvalues
%         of Df(0)
%         eigenvectors - corresponding eigenvectors
%         params - 2x1 vector in the form [mu, nu']
%         error - the error from the rigorous enclosure of the eigenpairs
% Output  coeff - multidimensional array of size (order+1) x (order+1) x 4
%         entry (i,j,k) corresponds to the coefficient k_{i-1,j-1} where in
%         the write up, k is either a,b,c,d and we subtract 1 from each
%         index to be consistent with the zero indexing in the text
%         Each matrix (:,:,i) will be upper triangular (upper left)
function coeff = calc_proj_coeff(order, eigenvalues, eigenvectors,params)
    % cast everything to an interval with appropriate radius
    coeff=zeros(order+1,order+1,4);
    e1=eigenvectors(:,1);
    e2=eigenvectors(:,2);
    lam1=eigenvalues(1);
    lam2=eigenvalues(2);
    
    Df0=JacSH(0,params(1),params(2));
    
    coeff(2,1,:)=e1;
    coeff(1,2,:)=e2;
    
% suborder=m+n and corresponds to the (m,n)th coefficient
    suborder=2;
    while suborder < order + 1 
        disp('Calculating terms of order:')
        disp(suborder)
        for j=suborder:-1:0
            i=suborder-j;
               
            Aij=(Df0-(i*lam1+j*lam2)*eye(4))^(-1);
            a=coeff(1:i+1,1:j+1,1);

            staraaij=starhat(a,a,i,j);
            staraaaij=tripstarhat(a,a,a,i,j);
            B=[0;0;0;-params(2)*staraaij+staraaaij];
            pij=Aij*B;
            coeff(i+1,j+1,:)=pij;
        end
  
        suborder=suborder+1;
    end
end


% gives the (mn)th element of the star hat cauchy product
function prod = starhat(a,b,m,n)
    sum=0;
    for i=0:m
        for j=0:n
            if (i == 0 &&  j == 0) || (i==m && j == n)
                % do nothing 
            else
                sum=sum+a(m-i+1,n-j+1)*b(i+1,j+1);
            end
        end
    end
    prod=sum;
end


% gives the (mn)th element of the star hat cauchy product between three
% elements
function prod = tripstarhat(a,b,c,m,n)
    sum=0;
    for j=0:m
        for k=0:n
            for l=0:j
                for q=0:k
                    if j == 0 && k == 0
                        % do nothing
                    elseif  l ==m && q == n 
                        % do nothing
                    elseif j==m && k ==n && l==0 && q == 0
                        % do nothing
                    else
                        sum=sum+a(m-j+1,n-k+1)*b(j-l+1,k-q+1)*c(l+1,q+1);
                    end
                end
            end
        end
    end
    prod=sum;
end

 



