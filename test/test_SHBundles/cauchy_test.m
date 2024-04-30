syms a [4 4] 
syms b [4 4] 


seq_a.coeff = a; 
seq_a.index_lb = -2; 
seq_a.dim = 2; 
seq_a.order = 1; 

seq_b.coeff = b; 
seq_b.index_lb = 0; 
seq_b.dim = 2; 
seq_b.order = 3; 

order = 1; 

test = Cauchy2(seq_a, seq_b, order);