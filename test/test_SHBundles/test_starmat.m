
v00 = [1;0]; 
v01 = [1;1]; 
v10 = [2;0];
v11 = [0;0];

v_coeff = zeros(2, 2, 2); 
v_coeff(:, 1, 1) = v00; 
v_coeff(:, 1, 2) = v01; 
v_coeff(:, 2, 1) = v10; 

A00 = [1,2;3,4]; 
A01 = [1,1;1,1]; 
A10 = [2,1;2,1]; 

A_coeff = zeros(2,2,2,2); 
A_coeff(:,:,1,1) = A00; 
A_coeff(:,:,1,2) = A01; 
A_coeff(:,:,2,1) = A10; 

Av00 = starhatMat(A_coeff,v_coeff,  0, 0);

Av10 = starhatMat(A_coeff, v_coeff, 1, 0);

Av01 = starhatMat(A_coeff, v_coeff, 0, 1);

convn(A_coeff, B_coeff)