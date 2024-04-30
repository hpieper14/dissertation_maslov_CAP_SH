% this function gets the value of the pulse (4 dim) at desiredTime 

function [t_0,ic] = getpulseIC(params, desiredTime, currentTime, currentIC) 
    options=odeset('MaxStep',0.001, 'RelTol', 1e-10);
    [t,P]=ode45(@(t,P)SH_comp(t,P,params),[currentTime desiredTime],currentIC, options);
    
    t_0 = desiredTime;
    ic = P(end, :);
end
