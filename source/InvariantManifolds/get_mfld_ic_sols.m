function solutions=get_mfld_ic_sols(num_sols,params,ics,end_time)
    solutions=cell(num_sols,4);
    
    for i = 1:num_sols
        tspan=[0 end_time];
        [st,sx]=ode45(@(st,sx) SH_comp(st,sx, params), tspan, ics(i,:));
        solutions(i,1)={st};
        solutions(i,2)={sx};
        solutions(i,3)={vecnorm(sx')'};
        
        odd_deriv=[sx(:,2),sx(:,4)];
        
        solutions(i,4)={vecnorm(odd_deriv.')};
    end



end