% This function performs runs the Taylor integrator to compute [phi, E^u_-]
% on the interval [orig_ic, orig_ic + stepSize] 
function [t,y] = stepTaylorIntegrator(orig_ic, orig_time, stepSize, params, order)
    time = orig_time:.01:orig_time+stepSize;
    coeffs =  frameTaylorIntegrator(orig_ic, order, params);
    % have to shift time around.
    tempTime = [0, stepSize];
    [t,y] = taylorCoeffToFunction(coeffs, tempTime);
    t = time; 
end
