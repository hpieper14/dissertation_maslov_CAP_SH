% PULSESOLUTION  Class object for the pulse solution.
classdef PulseSolution 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Class properties.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties( Access = public )
        
        % struct containing the fourier coefficients and associated parameters for Newtons method    
        fourier = [];         
        % parameters for the vector field
        vfParams = [];       
        % L for [-L,L] time interval for approximation
        time = [];         
        % struct containing normal form data and branch value
        normalForm = [];        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Constructor.                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    
    methods( Access = public, Static = false )
        
        function S = PulseSolution(fourier, vfParams, normalForm, time)
         
            % Store approximation order.
            S.fourier = fourier;
        
            % Store parameters for the vector field.
            S.vfParams = vfParams;

            % Store value of normal form branch 
            allowedBranchValues = [0,pi];
            if (ismember(normalForm.branch, allowedBranchValues) == 1) == 0
                error('The branch value is invalid. Please choose branch value 0 or pi.')
            end
            S.normalForm = normalForm;
            
            % Store temporal domain, should be a positive number for domain
            % [-L, L] 
            S.time = time;                     
            
        end
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Methods.                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods( Access = public, Static = false )
        % main function to obtain pulse approximation and unstable
        % eigenvalues associated to it
        S = mainPulse(S)
        
        % Normal form solution via Burke and Knobloch 
        S = BKNormalForm4d_halfline(S)

        % Get normal form solution with Dirichlet/Neumann BC 
        S = trimNFSol_halfline(S);

        % Get fourier coefficients from normal form solution 
        S = getFullFourierCoeffs(S, sol, interval_type)

        % Perform Newton's method 
        S = Newton_halfline(S); 
        S = Newton(S)

        % get pulse derivative initial conditions in symplectic coordinates
        S = getPulseDerivIC(S, t_0)
        
        % helper methods
        solution = getFunctionFromFourierCoeffs(S, coeffs, time_vec)
        jacobian = DFFourier(S, a)
        fourier_vf = fourierODE(S, a)    
                                    
    end   
    
end

