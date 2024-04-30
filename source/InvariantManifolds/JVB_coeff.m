function coeff_mat= JVB_coeff(m, tau)

mu = 0.2; 
nu = 1.6;


% Variable names. Here we interpret epsilon as a parameter (called xi).
% This is convenient for when we wish to perform continuation, but it is 
% not necessary, i.e., one does not have to include xi as an additional
% variable.   
variable_names = {'x1','x2','x3','x4'}; 

% Cell array which stores the components of the vectorfield. 
comp = cell(4,1); 
comp{1} = polynom(1,[0 1 0 0],variable_names);
comp{2} = polynom(1,[0 0 1 0],variable_names);
comp{3} = polynom(1,[0 0 0 1],variable_names); 
comp{4} = polynom([-2; -(mu+1); nu; -1],[0 0 1 0; 1 0 0 0; 2 0 0 0; 3 0 0 0],variable_names); 


% Construct vectorfield. 
g = Vectorfield(comp,[]);       % If we do not include xi as a parameter,
                                % then the syntax is g = Vectorfield(comp,[]) 
%g = g.reverse();                % We reverse time for convenience. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           II. Compute charts on the local (un)stable manifolds          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize a parameterization of the local stable manifold.   
eq_s = [0;0;0;0];                                 % Equilibrium (saddle-point).
n_s = 2;                                        % Dimension of the stable manifold.
K_s = [m m];                                  % Order of parameterization (initial "guess"). 
Q = Parameterization(g,n_s,K_s,eq_s,'stable');  % Initialize parameterization object.

% Refine the parameterization with Newton's method. 
Q = Q.rescale([tau tau]); % Need to rescale so that Newton is successful.  
Q = Q.compute_manifold(); 

coeff_mat=Q.coeffs;

% proj=[1,2,4];
% col='r';
% figure 
% Q.plot(proj,col)
end