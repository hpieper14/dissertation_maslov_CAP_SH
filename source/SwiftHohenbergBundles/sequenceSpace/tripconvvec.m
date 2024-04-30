function [vec] = tripconvvec(a,b,c,order, side)
    vec = [];
    if side == 1
        if (max(size(a)) == 3*order+1) == 0
            a = [a, zeros(1, 2*order)];
            b = [b, zeros(1, 2*order)];
            c = [c, zeros(1, 2*order)];
        end
        for k = 0:order 
            sum = 0;
            for k1 = -order:order
                for k2 = -order:order
                    k3 = k-k1-k2;
                    sum = sum+a(abs(k1)+1)*b(abs(k2)+1)*c(abs(k3)+1);
                end
            end
            vec = [vec, sum];
        end
    else
         if (max(size(a)) == 6*order+1) == 0
            a = [zeros(1,2*order), a, zeros(1, 2*order)];
            b = [zeros(1,2*order), b, zeros(1, 2*order)];
            c = [zeros(1,2*order), c, zeros(1, 2*order)];
         end
         for k = -order:order
             sum = 0;
             for k1 = -order:order
                 for k2 = -order:order
                     k3 = k - k1 - k2;
                     m = k1+3*order+1;
                     n = k2+3*order+1;
                     p = k3+3*order+1;
                     sum = sum+a(m)*b(n)*c(p);
                 end
            end
            vec = [vec, sum];
         end
    end
                     
        
    
end