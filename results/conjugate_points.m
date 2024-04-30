clear all 
close all 

% nonsnaking region
mu = 0.05; 
nu = 1.6; 
generate_results(0, mu, nu)
generate_results(pi, mu, nu)

% snaking region 
% TODO Add comment to readme about scaling for stable mu
mu = .2; 
nu = 1.6; 
generate_results(0, mu, nu)






function [S,C] = generate_results(branch, mu, nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normal form branch -- set to 0, pi
normalForm.branch = branch; 

% vector field parameters
vfParams.nu = nu;
vfParams.mu = mu;
vfParams.lambda = 0; 

% fourier approximation parameters
fourier.M = 1000; 
fourier.tol = 1e-14; 
fourier.order = 500; 
time = 150; 

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['We consider the pulse for parameter values nu=', ...
    num2str(vfParams.nu), ', mu=', ...
    num2str(vfParams.mu), ', and branch phi=', ...
    num2str(normalForm.branch), '.'])
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                PULSE SOLUTION APPROXIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

disp('Performing Newtons method to obtain a Fourier approximation of the pulse solution.')
disp(' ')

% initialize PulseSolution
S = PulseSolution(fourier, vfParams, normalForm, time);

% perform Newton's method
S = S.mainPulse();


% Plot normal form approximation and pulse approximation obtained by
% performing Newton's method 

time = S.normalForm.time; 
sol = S.normalForm.sol(:,1);
full_sol = S.getFunctionFromFourierCoeffs(S.fourier.full_coeff_from_half_newton, "full");
if S.vfParams.mu == .2 
    sol = 1/3.*sol; 
    full_sol(:,2) = 1/3.*full_sol(:, 2);
end

figure 
hold on 
plot(full_sol(:, 1), full_sol(:, 2), color = 'b', LineWidth=1.25)
plot([-flip(time), time], [flip(sol), sol], Color  = 'r', LineWidth=1.25)
%legend('Solution after Newtons method, $\bar \varphi^{(N)}$', 'Solution via Normal Form Eqn, $u_\phi$', Interpreter = 'latex')
legend('$\bar \varphi^{(N)}$', '$u_\phi$', Interpreter = 'latex', fontsize = 15)
title('Pulse Approximations, $\phi =$ ' + string(normalForm.branch) + ...
    ", $\mu =$ " + num2str(vfParams.mu) + ", $\nu =$ " + ...
    num2str(vfParams.nu), Interpreter = 'latex')
xlabel('$x$', Interpreter = 'latex', fontsize = 14)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 CONJUGATE POINT COMPUTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp(['Integrating the non-autonomous system (linearized about the pulse) ' ...
    'to obtain the basis solutions for $E^u_-(x, 0)$.'])
disp(' ')

% set parameters
conjPts.L = 60; 
C = ConjugatePoints(conjPts, vfParams);
C.Euminus.normalize = 1;
C.Euminus.refPlane = [1,4]; 

% compute conjugate points
[S,C] = C.mainConjPts(S);

%%%%%%%%%% Plot basis vectors for $E^u_-(x;0)$ %%%%%%%%%%%%%%%%%%%%%%%
time = C.Euminus.timeVec;
basis1 = reshape(C.Euminus.frame(:, 1,:), [4, max(size(time))]);
figure
tiledlayout(4,1)
nexttile
plot(time, basis1(1,:))
title('$W_1^{-,u}$, First Basis Vector for $E^u_-(x; 0)$, $\phi =$ ' + string(normalForm.branch) + ...
    ", $\mu =$ " + num2str(vfParams.mu) + ", $\nu =$ " + ...
    num2str(vfParams.nu), Interpreter = 'latex')
nexttile
plot(time, basis1(2,:))
nexttile
plot(time, basis1(3,:))
nexttile
plot(time, basis1(4,:))
xlabel('$x$', Interpreter = 'latex')

basis2 = reshape(C.Euminus.frame(:, 2,:), [4, max(size(time))]);
figure
tiledlayout(4,1)
nexttile
plot(time, basis2(1,:))
title('$W_2^{-,u}$, Second Basis Vector (Derivative of the Pulse) for $E^u_-(x; 0)$, $\phi =$ ' ...
    + string(normalForm.branch) + ...
    ", $\mu =$ " + num2str(vfParams.mu) + ", $\nu =$ " + ...
    num2str(vfParams.nu), Interpreter='latex')
nexttile
plot(time, basis2(2,:))
nexttile
plot(time, basis2(3,:))
nexttile
plot(time, basis2(4,:))
xlabel('$x$', Interpreter = 'latex')


phi1 = C.Euminus.phi1; 
phi2 = C.Euminus.phi2; 

figure
tiledlayout(1, 2)
nexttile 
plot(time, phi1(:,2))
title('Pulse solutions obtained via integrating, $\phi =$ ' ...
    + string(normalForm.branch) + ...
    ", $\mu =$ " + num2str(vfParams.mu) + ", $\nu =$ " + ...
    num2str(vfParams.nu), Interpreter='latex')
xlabel('$x$', Interpreter = "latex")
nexttile
plot(time, phi2(:,2))
xlabel('$x$', Interpreter = "latex")


%%%%%%%%% print eigenvalues and conjugate point locations %%%%%%%%% 
dets = C.conjPts.dets{:, 3}; 
time = C.conjPts.dets{:,2}; 

disp('Eigenvalue approximation via Fourier modes: ')
if size(S.fourier.unstable_eigs) == 0 
    disp('None found.')
else
    disp(S.fourier.unstable_eigs)
end

conj_ind = find([0, diff(sign(dets))] ~= 0);
conj_pts = time(conj_ind);

disp('Conjugate points locations (x): ')
if isempty(conj_pts)
    disp('None found.')
else
    disp(conj_pts)
end

disp("We found " + string(max(size(S.fourier.unstable_eigs))) + " eigenvalue(s) and "...
    + string(length(conj_pts)) + " conjugate points for the pulse with phase " ...
    +'phi = ' + string(normalForm.branch))
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% PLOT PULSE PROFILES, PLUCKER EMBEDDING, DET PLOTS %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Deterimant plot
figure 
plot(time, dets, LineWidth=1.25)
title('Determinant $\det(A)$, $\phi =$ ' + string(normalForm.branch) + ...
    ", $\mu =$ " + num2str(vfParams.mu) + ", $\nu =$ " + ...
    num2str(vfParams.nu), Interpreter='latex')
xlabel('$x$', Interpreter = 'latex', fontsize = 14)
ylabel('$\det A(x)$', Interpreter = 'latex', fontsize = 14)
if  mu == .2 
    axis([-conjPts.L, conjPts.L, 0, .8])
end

if vfParams.mu == .2
    c = "#0072BD";
elseif normalForm.branch == pi 
    c = "#77AC30"; 
else
    c = "#A2142F"; 
end

% pulse profile figures first shown in introduction
time = S.normalForm.time; 
full_sol = S.getFunctionFromFourierCoeffs( ...
    S.fourier.full_coeff_from_half_newton, "full");
if S.vfParams.mu == .2 
    full_sol(:,2) = 1/3*full_sol(:,2);
end
figure 
plot(full_sol(:, 1), full_sol(:, 2), color = c, LineWidth=1.25)
xlabel('$x$', Interpreter = 'latex', fontsize = 14)

frame_trajectory = C.Euminus.frame;
time_length = length(frame_trajectory(1,1,:));
plucker_trajectory = zeros(time_length,6);
for i=1:time_length
    v1 = frame_trajectory(:,1,i);
    v2 = frame_trajectory(:,2,i);
    plucker_trajectory(i,:)= plucker_coords(v1,v2);
end

plot_time = (max([1,ceil(time_length*.0)]):time_length);
coord_12 = plucker_trajectory(plot_time,1);
coord_13 = plucker_trajectory(plot_time,2);
coord_14 = plucker_trajectory(plot_time,3);

% coordinates for sandwich plane train
[x,y] = meshgrid(-1:0.1:1); 
z = zeros(size(x, 1));

if vfParams.mu == .2
    trajectory_color = .7*[0 0.4470 0.7410];
elseif normalForm.branch == pi 
    trajectory_color = .7*[0.4660 0.6740 0.1880];
else 
    trajectory_color = .7*[0.6350 0.0780 0.1840];
end

% Plot of trajectory and sandwich plane train in Plucker coordinates
set(gcf, 'Renderer', 'Painters');
figure
grid on
hold  on 
plot3(coord_12,coord_13,coord_14, LineWidth=1.5, color = trajectory_color);
surf(x, y, z, FaceAlpha = 0.15, EdgeColor="none", FaceColor = "#7E2F8E");
xlabel('$P_{12}$', Interpreter = 'latex', FontSize=14)
ylabel('$P_{13}$', Interpreter = 'latex', FontSize=14)
zlabel('$P_{14}$', Interpreter = 'latex', FontSize=14) 
view(-20, 10)

end

function [coord] = plucker_coords(v,w)
%UNTITLED2 Computes the wedge product of two 4-vectors, then normalizes
   % Coordinates given by 
    % 1 --  x_1 /\ x2 
    % 2 --  x_1 /\ x3 
    % 3 --  x_1 /\ x4 
    % 4 --  x_2 /\ x3 
    % 5 --  x_2 /\ x4
    % 6 --  x_3 /\ x4 
    
    coord = 0*(1:6);
    wedge_index = zeros(6,2);
    wedge_index(1,:) =[1 2];
    wedge_index(2,:) =[1 3];
    wedge_index(3,:) =[1 4];
    wedge_index(4,:) =[2 3];
    wedge_index(5,:) =[2 4];
    wedge_index(6,:) =[3 4];
    
    for i = 1:6
        index_right = wedge_index(i,1);
        index_left = wedge_index(i,2);
    
        coord(i) = v( index_left )*w(index_right ) - v( index_right )*w(index_left );
    end
    coord = coord/norm(coord,2);

end

 