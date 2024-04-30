max_order = 1; 
lam = [1,2]; 

coeff = lam_dot_alpha_mat(lam, max_order); 

by_hand = [-6, -4, -2, 0; 
            -5, -3, -1, 1;
            -4, -2, 0, 2; 
            -3, -1, 1, 3]; 

assert(isequal(coeff, by_hand))


lam = [1 + 1i*2, 1 - 1i*2];

coeff = lam_dot_alpha_mat(lam, max_order);

by_hand = [-4, -3-2*1i, -2-4*1i, -1-6*1i; 
           -3+2*1i, -2, -1-2*1i, -4*1i; 
           -2+4*1i, -1+2*1i, 0, 1-2*1i; 
           -1+6*1i, 4*1i, 1+2*1i, 2];

assert(isequal(coeff, by_hand))

