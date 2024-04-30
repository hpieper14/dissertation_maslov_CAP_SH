function Aa_norms = get_Aa_norms(nu, x, params, mflds)
    m = params.cheb.order;

    nu_power=nu.^(0:m-1);

    a1_nu_norm=sum(abs(x.a1)'.*nu_power');
    a2_nu_norm=sum(abs(x.a2)'.*nu_power');
    a3_nu_norm=sum(abs(x.a3)'.*nu_power');
    a4_nu_norm=sum(abs(x.a4)'.*nu_power');
    
    Aa_norms.a1 = a1_nu_norm;
    Aa_norms.a2 = a2_nu_norm;
    Aa_norms.a3 = a3_nu_norm;
    Aa_norms.a4 = a4_nu_norm;
    
    disp(['The nu-norm of a1 is = ',num2str(a1_nu_norm)])
    disp(['The nu-norm of a2 is = ',num2str(a2_nu_norm)])
    disp(['The nu-norm of a3 is = ',num2str(a3_nu_norm)])
    disp(['The nu-norm of a4 is = ',num2str(a4_nu_norm)])

    Am = (DF_homoclinic(x,params,mflds))^(-1);
    
    A_norms=zeros(7,7);
    A_norms(1:3,1:3)=abs(Am(1:3,1:3));
    for j=4:7
        A_norms(1:3,j)=max(abs(Am(1:3,3+1+(j-4)*m:3+(j-3)*m).*repmat(nu.^-(0:m-1),3,1)),[],2);
    end
    for i=4:7
        A_norms(i,1:3)=sum(abs(Am(3+(i-4)*m+1:3+(i-3)*m,1:3)).*repmat(nu_power',1,3),1);
    end
    
    K_A=zeros(7,7);
    for i=4:7
        for j=4:7
            K_A(i,j)=max(mag(sum(abs(Am(3+(i-4)*m+1:3+(i-3)*m,3+(j-4)*m+1:3+(j-3)*m)).*repmat(nu_power,m,1)').*(nu.^-(0:m-1))));

        end
    end
    A_norms(4:7,4:7)=max(K_A(4:7,4:7),sup(1/(2*m)));
    disp('The norms of the sub operators of A are ')
    disp(A_norms)
    
    Aa_norms.A = A_norms;
end
