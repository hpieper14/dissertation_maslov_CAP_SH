% need to fix order so that the returned vector has length order -- need to
% then adjust the use in the proof of the homoclinic orbit existence
% the coefficients of order N are stored in the N+1 entry
function vec = chebstar2(a,b,order)
    N=max(size(a));
    M=max(size(b));
    
    if (max(size(a)) == order) == 0
        a=[a,zeros(1,order-N)];
    end
    if (max(size(b)) == order) == 0
        b=[b,zeros(1,order-M)];
    end
    
        
    vec=[];
    for k=0:order-1
        dubsum=0;
        for i=-order+1:order-1
            l=k-i;
            if abs(l)>order-1
                % do nothing 
            else
                dubsum=dubsum+a(abs(i)+1)*b(abs(l)+1);
            end
        end
        vec=[vec,dubsum];
    end

end
