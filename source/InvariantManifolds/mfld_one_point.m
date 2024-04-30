function val = mfld_one_point(phi1, phi2, coeff, params)
    order=params.mfld.order;
    val=zeros(4,1);
    for n=0:order
        for m=0:n
            point=reshape(coeff(n-m+1,m+1,:),[4,1]);
            val=val+point.*(phi1+1i*phi2)^(n-m)*(phi1-1i*phi2)^m;

        end
    end
       if norm(val-real(val)) > 1e-10
            msg = 'Error occurred. The manifold is complex valued.';
            error(msg);
       end
end 

