function injective = check_A_injective(x,params,mflds)
    m = params.cheb.order;

    DF=DF_homoclinic(x,params,mflds);
    Am=DF^(-1);

    disp('Norm of I - ADF: ')
    disp(norm(mag(eye(4*m+3)-Am*DF)))

    if norm(mag(eye(4*m+3)-Am*DF))>=1
        disp('Stop! The matrix Am is not injective.');
        injective = 0;
        return
    else
        disp('Good to go! The matrix Am is injective.');
        injective=1;
    end
end
