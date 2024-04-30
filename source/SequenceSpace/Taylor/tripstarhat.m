% we format a, b as N+1 square matrices
% gives the (mn)th element of the star hat cauchy product
function prod = tripstarhat(a,b,c,m,n)
    sum=0;
    for j=0:m
        for k=0:n
            for l=0:j
                for q=0:k
                    if j == 0 && k == 0
                        % do nothing 
                    elseif j==m && l ==m && q == n && k == n
                        % do nothing 
                    elseif j==m && k ==n && l==0 && q == 0
                        % do nothing 
                    else
                        sum=sum+a(m-j+1,n-k+1)*b(j-l+1,k-q+1)*c(l+1,q+1);
                    end
                end
            end
        end
    end
    prod=sum;
end
