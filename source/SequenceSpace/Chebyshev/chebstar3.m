function vec = chebstar3(a,b,c,order)
    vec=[];
    for k=0:order-1
        %disp('-------------------------------')
        %disp(['order',num2str(k)])
        tripsum=0;
        for i=-order+1:order-1
            for j=-order+1:order-1
                l=k-j-i;
                if abs(l) > order-1
                    % do nothing 
                else 
                    tripsum=tripsum+a(abs(i)+1)*b(abs(j)+1)*c(abs(l)+1);
         %           disp(['Signed indices: ', num2str(i),' ', num2str(l),' ' , num2str(j)])
          %          disp(['order=? ',num2str(i+l+j)])
                end
                
            end
        end
        vec=[vec,tripsum];
    end 

end