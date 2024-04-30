% we format a, b as N+1 square matrices
% gives the (mn)th element of the star hat cauchy product
function prod = starhat(a,b,m,n)
    % if m+n > size(a)
    %     error('Desired index exceeds truncation')
    % end
    sum=0;
    for i =0:m
        for j=0:n
            if i == 0 &&  j == 0
                sum=sum+0;
            elseif i == m && j == n
                sum=sum+0;
            else
            sum=sum+a(m-i+1,n-j+1)*b(i+1,j+1);
            end
        end
    end
    prod=sum;
end