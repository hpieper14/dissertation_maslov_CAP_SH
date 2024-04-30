function sum = dubstarhat_1d(a,b,c,m)
    sum = 0; 
    for i = 0:m 
        for j = 0:m-i
            k = m - i - j;
            if i == m || j == m 
                % do nothing
            elseif i == 0 || j == 0 
                % do nothing
            else 
                sum = sum + a(i+1)*b(j+1)*c(k+1);
            end
        end
    end
end