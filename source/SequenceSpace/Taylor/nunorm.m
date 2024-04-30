% computes the nu norm. In this program, N=1
% Inputs: vec - object to obtain norm of  
%         nu - weight
%         N - should be 1  
% Outputs norm - nonzero scalar
function norm = nunorm(vec,N,nu)
    if nu < 0 
        disp('Please choose a positive weight')
    end
    norm= (0);
    if N==1
        if min(size(vec))==1==false 
            disp('N and object size do not agree')
        end
        s=max(size(vec));
        for i=0:s-1
            norm=norm+abs(vec(i+1))*nu^i;
        end
    else 
        A =  (zeros(1,N)); 
        for j=1:N
            norm= (0);
            s=max(size(vec));
            for i=0:s-1
                norm=norm+abs(vec(i+1))*nu^i;
            end
            A(j)=norm;
        end
    norm = max(A);
    end
end