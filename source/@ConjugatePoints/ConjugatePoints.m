% CONJUGATEPOINTS ConjugatePoints Class. Used to compute the conjugate
% points associated to PulseSolution object S 

classdef ConjugatePoints 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Class properties.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties( Access = public )

        conjPts = [];         % struct containing data associated to the conjugate points  
        vfParams = [];         % Parameters for the vector field
        Euminus = []; 
        Esplus = [];
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Constructor.                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    
    methods( Access = public, Static = false )
        
        function C = ConjugatePoints(conjPts, vfParams)
         
            % Store approximation order.
            
            % Store parameters for the vector field.
            C.vfParams = vfParams;
            C.conjPts = conjPts; 
                                    
            
        end
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Methods.                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods( Access = public, Static = false )

        [S, C] = mainConjPts(C, s)
        
        mat = Binf(C); 
        [vectors, values]= getBinfEigs(C);
        isLagrangian = verifyLagrangian(C, frame);

        C = generateEuFrame(C);
            
        vals= calculateDeterminant(C);

        C = getConjPtLocations(C, vals); 

        [frame, time]= makeFrame(C, basis_1, basis_2)
        C = get4DimIC(C, S)




                                    
    end   
    methods( Access = private, Static = false )

        f1 = nonautonODE(C, T, Y)
        f1 = nonautonODENormalized(C, T, Y)
        f1 = stableBasisNonautonODE(C,T,Y, eval)
                
        
                                    
    end   
    
end

