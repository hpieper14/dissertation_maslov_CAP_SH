% GENERATEEUFRAME  Generates the frame for $E^u_-(x;0)$ for various x
% values and saves it to the ConjugatePoints object C 
% 
%   C = generateEuFrame(C)
%   C = C.generateEuFrame()
function C = generateEuFrame(C)
    L = C.conjPts.L; 
    orig_ic = C.Euminus.pulse_ic';
    deriv_ic = C.Euminus.pulse_deriv_ic';
    
    % inialize initial condition for the basis solution of $E^u_-$ that
    % is not constructed from the derivative of the pulse
    [vectors, values]= getBinfEigs(C);
    v1=real(vectors.u(:,1));    

    if C.Euminus.normalize == 1 
        deriv_ic = deriv_ic./vecnorm(deriv_ic);
        v1 = v1./vecnorm(v1);
    end

    new_ic1 = [orig_ic; v1];
    new_ic2 = [orig_ic; deriv_ic]; 

    [full_phi1, basis_1] = C.integrateDE(new_ic1,-L,L);
    [full_phi2, basis_2] = C.integrateDE(new_ic2,-L,L);

    C.Euminus.phi1 = full_phi1; 
    C.Euminus.phi2 = full_phi2; 
    
    [C.Euminus.frame, C.Euminus.timeVec] = C.makeFrame(basis_1, basis_2);
   
end    