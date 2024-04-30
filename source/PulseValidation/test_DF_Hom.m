function val=test_DF_Hom(epsilon)

    [x,params, mflds]=init_sol_val_tests();
    

    approx_df = approx_DFHom(x, epsilon,params,mflds);
    DF = DF_homoclinic(x, params, mflds);
    
    val=approx_df-DF;

    if max(val,[],'all')>1e-4
        disp('There might be something wrong with DF.')
    end

end
