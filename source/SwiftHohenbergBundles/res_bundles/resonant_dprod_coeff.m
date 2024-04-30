% this assumes growth rate corresponds to the first eigenvalue in the
% vector eigenvalues.stable 
function [w_coeff] = resonant_dprod_coeff( ...
    eigenvalues, v_coeff, order, growth_rate)
    offset = 3;

    % assert order < max_order(v_coeff) - 2 (I think) 

    v0_1 = v_coeff(0+offset, 0 + offset); 
    lam = eigenvalues.s; 

    vI = dprod_coeff_v_I(v_coeff, order, lam, growth_rate); 
    vII = dprod_coeff_v_II(v_coeff, order, lam, growth_rate); 
    vIII = dprod_coeff_v_III(v_coeff, order, lam, growth_rate); 

    w_coeff = zeros(order + offset, order + offset); 

    % will need to change this assignment based on the chosen growth rate
    w_coeff(-2 + offset, 0 + offset) = 1; 

    lama_w = zeros(order + offset, order + offset); 
    lama2_w = zeros(order + offset, order + offset); 
    lama3_w = zeros(order + offset, order + offset); 
    
    i = -2; 
    j = 0; 
    lama_w(i + offset, j + offset) = dot(lam, [i,j]);
    lama2_w(i + offset, j + offset) = dot(lam, [i,j])^2;
    lama3_w(i + offset, j + offset) = dot(lam, [i,j])^3;

 
    for alpha = -1:1:order 
        for a1 = -2:1:alpha 
            a2 = alpha - a1; 
            if a2 > order
                continue 
            else 
                  
                a_vec = [a1, a2];
                d0 = growth_rate^3 - growth_rate; 
                d1 = dot(a_vec, lam)*(6*growth_rate^2 - 2);
                d2 = growth_rate*dot(a_vec, lam)^2; 
                d3 = dot(a_vec, lam)^3; 
    
                frac = -(v0_1*(d0 + d1 + d2 + d3))^(-1); 
    
                term1 = starhat(w_coeff, vI, a1, a2); 
                term2 = starhat(lama_w, vII, a1, a2); 
                term3 = starhat(lama2_w, vIII, a1, a2); 
                term4 = starhat(lama3_w, v_coeff, a1, a2); 
                
                wa1a2 = frac*(term1 + term2 + term3 + term4);
    
                w_coeff(a1 + offset, a2 + offset) = wa1a2;
    
                lama_w(a1 + offset, a2 + offset) = dot(lam, [a1, a2])*wa1a2; 
                lama2_w(a1 + offset, a2 + offset) = dot(lam, [a1, a2])^2*wa1a2; 
                lama3_w(a1 + offset, a2 + offset) = dot(lam, [a1, a2])^3*wa1a2; 
            end
        end
    end
end
