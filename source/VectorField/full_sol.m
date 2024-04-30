function [full_phi, full_v] = full_sol(u_ic, s_ic,L1, L2, params, vec)
    % integrate to obtain left half of solution
    Lic=[u_ic, vec];
    [Lphisol,Lvsol]=integrateDE(Lic, -L1,0, params);
    % integrate to obtain right half of solution
    Ric=[s_ic, vec];
    [Rphisol,Rvsol]=integrateDE(Ric, L1,0, params);
    
    if (L2==0)==0
        [Lphitail,Lvtail]=integrateDE(Lic, -L1, -L2, params);
        [Rphitail,Rvtail]=integrateDE(Ric, L1, L2,params);
        full_phi=[flip(Lphitail);Lphisol;flip(Rphisol); Rphitail];
        full_v=[flip(Lvtail);Lvsol;flip(Rvsol);Rvtail];
    else
        full_phi=[Lphisol;flip(Rphisol)];
        full_v=[Lvsol;flip(Rvsol)];
    end
end