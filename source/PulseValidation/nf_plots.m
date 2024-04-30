clear all 
clear figure 

params=[.05,1.6]; %[mu, nu]

Df0=JacSH(0,params(1),params(2));
[V,D]=eigs(Df0);

tau=1/4; 
Vscale=tau*V;

stabeigs=[D(1,1),D(2,2)];
stabvec=Vscale(:,1:2);

uneigs=[D(3,3),D(4,4)];
unvec=Vscale(:,3:4);


order=20;



% unstable
uncoeff=calc_proj_coeff(order,uneigs,unvec,params);

% stable 
stcoeff=calc_proj_coeff(order,stabeigs,stabvec,params);

phi=pi/2;

x=-20:1/10:20;
y=sh_nf_burke_knobloch(x,phi,params);

 
time=10;

ics=sh_nf_burke_knobloch(0,phi,params);

figure
hold on
plot_manifold(uncoeff,order, 'r');
plot_manifold(stcoeff,order, 'g');
plotsoln(params, ics,time);
plot3(y(1,:),y(2,:),y(4,:));
title('Manifolds and NF Solution')
hold off 

figure 
acloserlook(params, sh_nf_burke_knobloch(0,phi,params), time);



q_2=3/4-19*params(2)^2/18;
x=-55:1/10:55;
y = sqrt(-2.*params(1)/q_2).*sech(x.*sqrt(params(1))/2).*cos(x+phi);

figure 
plot(x,y);
title('Plot of the NF solution')

y=sh_nf_burke_knobloch(x,phi,params);
energies=get_energy_levels(params, y');
N=max(size(energies));
x=linspace(-40,40,N);

figure 
plot(x,energies)
title('Energy Levels of NF Solution')


% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val=get_energy_levels(params, soln)
    N = max(size(soln));
    integrand=zeros(1,N);
    
    for i=1:N
        u=soln(i,1);
        u_x=soln(i,2);
        u_xx=soln(i,3);
        mu=params(1);
        nu=params(2);
        integrand(i)=mu/2*u^2-nu/3*u^3+1/4*u^4 + 1/2*((u_xx)^2-2*(u_x)^2+u^2);
    end
    energy_vec=cumtrapz(integrand);
    
    u=soln(1,1);
    u_x=soln(1,2);
    u_xx=soln(1,3);
    E=mu/2*u^2-nu/3*u^3+1/4*u^4 + 1/2*((u_xx)^2-2*(u_x)^2+u^2);
    
    Guess=E*norm(soln(1,:));
    val=[Guess,energy_vec];
end

function plots=acloserlook(params,ic,time) 
    tspan=[-time time];
    [t, x] = ode45(@(t,x) SH_comp(t,x, params), tspan, ic);
    disp(size(x))
    N = max(size(t));
    norms=zeros(1,N);
    for i=1:N 
        norms(i)=norm(x(i,:));
    end
    
    tiledlayout(3,1)
    nexttile
    plot3(x(:,1), x(:,2), x(:,4));
    title('Solution in 3-dim space')
    nexttile
    plot(t, norms);
    title('Norm of Soln over time')
    nexttile
    plot(t,x(:,1));
    title('Trajectory of x_1')
    plots=0;
end 


function soln = plotsoln(params, ic, time)
    tspan=[-time time];
    [t, x] = ode45(@(t,x) SH_comp(t,x, params), tspan, ic); 
    soln=plot3(x(:,1), x(:,2), x(:,4));
end 


    
 


function plots=plot_coeff(coeff,order)
    normorder=zeros(order+1,order+1);
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
                normorder(i+1,j+1)=normpoint;
                point=[suborder;normpoint];
                plotpoints=[plotpoints,point];
            end
        end
        suborder=suborder+1;
    end
    
  plot(plotpoints(1,:), log(plotpoints(2,:)),'o');
  plots=0;
end
% params in the form (mu, nu)

% coeff mnth coeff is stored in p(m,n,:)
function coeff = calc_proj_coeff(order, eigenvalues, eigenvectors,params)
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
        for i=0:suborder
            for j=0:suborder-i
                if (i==0 && j==0) || (i == 0 && j == 1) || (i==1 && j == 0)
                  % do nothing 
                else 

                    Aij=(Df0-(i*lam1+j*lam2)*eye(4))^(-1);
                    a=coeff(1:i+1,1:j+1,1);

                    staraaij=starhat(a,a,i,j);
                    staraaaij=tripstarhat(a,a,a,i,j);
                    B=[0;0;0;-params(2)*staraaij+staraaaij];
                    pij=Aij*B;
                    coeff(i+1,j+1,:)=pij;
                end
            end
        end
        suborder=suborder+1;
    end
end



% we format a, b as N+1 square matrices
% gives the (mn)th element of the star hat cauchy product
function prod = starhat(a,b,m,n)
    sum=0;
    for i=0:m
        for j=0:n
            if i == 0 &&  j == 0
                sum=sum+0;
            elseif i == m && j == n
                sum=sum+0;
            else
                
                sum=sum+a(m-i+1,n-j+1)*b(i+1,j+1);
            end
        end
    end
    prod=sum;
end

% we format a, b as N+1 square matrices
% gives the (mn)th element of the star hat cauchy product
function prod = tripstarhat(a,b,c,m,n)
    sum=0;
    for j=0:m
        for k=0:n
            for l=0:j
                for q=0:k
                    if j == 0 && k == 0
                        sum=sum+0;
                    elseif j==m && l ==m && q == n && k == n
                        sum=sum+0;
                    elseif j==m && k ==n && l==0 && q == 0
                        sum=sum+0;
                    else
                        sum=sum+a(m-j+1,n-k+1)*b(j-l+1,k-q+1)*c(l+1,q+1);
                    end
                end
            end
        end
    end
    prod=sum;
end


