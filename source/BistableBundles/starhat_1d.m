function sum = starhat_1d(a,b,m)
    sum=0;
    for i=1:m-1
        sum = sum + a(m-i+1)*b(i+1);
    end
end