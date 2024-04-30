function coeffs = DFQbundle(params, mflds)
   Q = mflds.coeff.uu;
   
   Q1 = Q(:,1);
    
   order = params.mfld.order;
   q1q1 = star_1d(Q1, Q1, order);
   
   coeffs = zeros(2,2,order+1);
   coeffs(:,:,1) = [0,1; 1, 0];
    
    for i = 1:order 
            coeff21 =  3*params.b*(-Q1(i+1) + q1q1(i+1));
            coeffs(:,:,i+1) = [0,0; coeff21, 0];
    end
end