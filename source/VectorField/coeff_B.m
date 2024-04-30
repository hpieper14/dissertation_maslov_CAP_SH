function matrix = coeff_B(u1,lambda, mu, nu)
    fprimeu=-mu+2*nu*u1-3*u1^2;
    matrix = [0,0,0,1;0,0,1,-2;lambda-1+fprimeu,0,0,0;0,1,0,0];
end