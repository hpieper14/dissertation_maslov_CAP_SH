%diary Newton_chebyshev_BK

% Set initial seed for Newton's method 
[x,params, mflds] = init_sol_val_tests();

params.normalizeBasis = 0;

disp('Choice of Parameters')
disp(params)
disp(['Chebyshev order: ', num2str(params.cheb.order)])
disp(['Manifold order: ', num2str(params.mfld.order)])

disp(['Initial choice for phi1, phi2, psi: ', num2str(x.phi1), ', ', num2str(x.phi2), ', ', num2str(x.psi)]);


u_ic = [-0.190343484893489 ; -0.128234936941383 ; 0.193875175006471 ; 0.168311653800752];
s_ic = [-0.190343484893489 ; 0.128234936941383 ; 0.193875175006471 ; -0.168311653800752];

L1=params.L;


[V, D]=asym_un_vecs(params);
v1=real(V(:,1));

[full_phi1, basis_1]= full_sol(u_ic, s_ic, L1,0, params, v1);

% %%% Option to use shooting solution
%  sol=full_phi1(:,2:5);
%  time_vec=full_phi1(:,1);

%%% Option to use NF solution
sol=BK_nf_4dim(params,0,L1, 0);
time_vec=sol(:,1);
sol=sol(:,2:5);

 figure
 tiledlayout(4,1)
 nexttile
 plot(time_vec,sol(:,1))
 nexttile
 plot(time_vec, sol(:,2))
 nexttile
 plot(time_vec,sol(:,3))
 nexttile
 plot(time_vec,sol(:,4))
 title('Solution obtained via BK.')

cheb_coeff = get_cheb_coeffs(sol, params);
x.a1=cheb_coeff(:,1)';
x.a2=cheb_coeff(:,2)';
x.a3=cheb_coeff(:,3)';
x.a4=cheb_coeff(:,4)';

yo1 = chebcoeff_to_function(x.a1);
yo2 = chebcoeff_to_function(x.a2);
yo3 = chebcoeff_to_function(x.a3);
yo4 = chebcoeff_to_function(x.a4);

nozeros=abs(nonzeros(x.a1));
max_ord=my_order(max(nozeros));
min_ord=my_order(min(nozeros));

disp(['Decay in Chebyshev coefficients for first component: ', num2str(max_ord-min_ord)])

G=FHomoclinic(x, mflds, params);

% Perform Newton's method 

new_x = refine_cheb_orbit(x, mflds, params);


% Compute Bounds for Radii Polynomial 

verif = verify_homoclinic_orbit(params, mflds, new_x, 1.05);


%%%%%%%%%%%%%%%%%
% Generate Plots 
%%%%%%%%%%%%%%%%

y1 = chebcoeff_to_function(new_x.a1);
y2 = chebcoeff_to_function(new_x.a2);
y3 = chebcoeff_to_function(new_x.a3);
y4 = chebcoeff_to_function(new_x.a4);

dom = -1:.05:1;

 figure
 tiledlayout(4,1)
 nexttile
 plot(dom,y1)
 nexttile
 plot(dom,y2)
 nexttile
 plot(dom,y3)
 nexttile
 plot(dom,y4)
 title('Solution obtained via Newtons method.')


figure 
hold on 
plot(dom, y1, linewidth = 1.5, color = "#A2142F")
plot(dom(1), y1(1), 'o', 'MarkerFaceColor', "#A2142F")
plot(dom(end), y1(end), 'o', 'MarkerFaceColor', "#A2142F")
hold off
xlabel('$t$', Interpreter = 'latex', FontSize=14)
ylabel('$\varphi(t)$', Interpreter = 'latex', FontSize=14)

 
 
 dom = params.L*dom;
 
 figure 
 tiledlayout(4,1)
 nexttile
 hold on
 plot(time_vec,sol(:,1))
 plot(dom,yo1)
 plot(dom,y1)
 legend('BK Solution','Chebyshev Rep.', 'Refined Chebyshev Rep.')
 hold off
 nexttile
 hold on
 plot(time_vec,sol(:,2))
 plot(dom,yo2)
 plot(dom,y2)
 hold off
 nexttile
 hold on
 plot(time_vec,sol(:,3))
 plot(dom,yo3)
 plot(dom,y3)
 hold off
 nexttile
 hold on 
 plot(time_vec,sol(:,4))
 plot(dom,yo4)
 plot(dom,y4)
 hold off 
 
 
figure
tiledlayout(2,1)
nexttile
plot_coeff(mflds.unstable.coeffs, params.mfld.order);
title('Unstable Mfld Coeff.')
nexttile
plot_coeff(mflds.stable.coeffs, params.mfld.order);
title('Stable Mfld Coeff.')


mflds.pts.s=mfld_points(mflds.stable.coeffs, params);
mflds.pts.u=mfld_points(mflds.unstable.coeffs, params);

figure
hold on
surf(mflds.pts.s(:,:,1),mflds.pts.s(:,:,2),mflds.pts.s(:,:,4), 'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none');  
xlabel('x1');
ylabel('x2');
zlabel('x4');
surf(mflds.pts.u(:,:,1),mflds.pts.u(:,:,2),mflds.pts.u(:,:,4), 'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none');  
plot3(y1,y2,y4);
title('Stable and Unstable Manifolds');
legend('Stable','Unstable','Hom. Orbit');
hold off



figure
grid on 
hold on
colormap(spring)
surf(mflds.pts.s(:,:,1),mflds.pts.s(:,:,2),mflds.pts.s(:,:,4), mflds.pts.s(:,:,3), 'FaceAlpha',0.4, 'EdgeColor', '#554b1c');
xlabel('$x_1$', Interpreter = 'latex');
ylabel('$x_2$', Interpreter = 'latex');
zlabel('$x_4$', Interpreter = 'latex')
surf(mflds.pts.u(:,:,1),mflds.pts.u(:,:,2),mflds.pts.u(:,:,4),'FaceAlpha',0.4);  
plot3(y1,y2,y4, lineWidth = 2, Color="#A2142F");
plot3(y1(1), y2(1), y4(1), 'o', 'MarkerFaceColor', "#A2142F")
plot3(y1(end), y2(end), y4(end), 'o', 'MarkerFaceColor', "#A2142F")
%title('Validated Manifolds and $\varphi(x)$ Trajectory from ');
legend('Stable manifold','Unstable manifold','$\varphi(x)$', Interpreter = 'latex');
hold off



figure
grid on 
hold on
surf(mflds.pts.s(:,:,1),mflds.pts.s(:,:,2),mflds.pts.s(:,:,4), 'FaceAlpha',0.4, 'FaceColor', '#554b1c', edgeColor = "none");
xlabel('$x_1$', Interpreter = 'latex');
ylabel('$x_2$', Interpreter = 'latex');
zlabel('$x_4$', Interpreter = 'latex')
surf(mflds.pts.u(:,:,1),mflds.pts.u(:,:,2),mflds.pts.u(:,:,4),'FaceAlpha',0.4, 'FaceColor',"#7E2F8E", edgeColor = "none"); 
plot3(y1,y2,y4, lineWidth = 2, Color="#A2142F");
plot3(y1(1), y2(1), y4(1), 'o', 'MarkerFaceColor', "#A2142F")
plot3(y1(end), y2(end), y4(end), 'o', 'MarkerFaceColor', "#A2142F")
%title('Validated Manifolds and $\varphi(x)$ Trajectory from ');
legend('Stable manifold','Unstable manifold','$\varphi(x)$', Interpreter = 'latex');
hold off
