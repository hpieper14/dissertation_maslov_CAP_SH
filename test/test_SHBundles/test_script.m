
order = 2; 
a = []; 
a.index_lb = -1; 
a.coeff = [1,2,3; 1,1, 1; 4,5,6]; 
a.dim = 2; 
a.order = order; 

b = []; 
b.index_lb = -1; 
b.coeff = [4,1,2; 1,1,1; 1,2,3]; 
b.dim = 2; 
b.order = order;

test = cauchy(a, b, order); 
test2 = Cauchy2(a,b,order);
% computed by hand
astarb = [4,9,16; 5, 8, 0; 18, 0, 0];

assert(isequal(astarb, test2))    
assert(isequal(test2,test.coeff))
    


function test_conv_cauchy_0_indexing() 
    order = 2; 
    a = []; 
    a.index_lb = 0; 
    a.coeff = [1,2,3; 1,1, 1; 4,5,6]; 
    a.dim = 2; 
    a.order = order; 
    
    b = []; 
    b.index_lb = 0; 
    b.coeff = [4,1,2; 1,1,1; 1,2,3]; 
    b.dim = 2; 
    b.order = order;
    
    test = cauchy(a, b, order); 
    test2 = Cauchy2(a,b,order);
    % computed by hand
    astarb = [4,9,16; 5, 8, 0; 18, 0, 0];

    assert(isequal(astarb, test2))    
    assert(isequal(test2,test.coeff))
end 


function test_syntax_indices_in_bounds()
    offset = 3;
    order = 10; 
    
    v_coeff = rand(order + offset, order + offset);
    v_coeff(1,1) = 0; 
    v_coeff(1,2) = 0; 
    v_coeff(2,1) = 0; 
    v_coeff(2,2) = 0; 
    
    
    eigenvalues.stable = rand(1,2); 
    eigenvalues.unstable = rand(1,2);
    
    growth_rate = eigenvalues.unstable(1); 
    
    
    order_to_compute = 8;
    [coeff] = resonant_dprod_coeff( ...
        eigenvalues, v_coeff, order_to_compute, growth_rate);
end
