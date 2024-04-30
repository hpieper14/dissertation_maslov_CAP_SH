function dydt=SH_comp(t,x,params)
    dydt= zeros(4,1);

    dydt(1)=x(2);
    dydt(2)=x(3);
    dydt(3)=x(4);
    dydt(4)= -2.*x(3) - (params.mu+1).*x(1) + params.nu.*x(1).^2 - x(1).^3;
end