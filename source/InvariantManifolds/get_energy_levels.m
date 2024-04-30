function val=get_energy_levels(params, soln)
    N = max(size(soln));
    energy_vec=zeros(1,N);

    for i=1:N
        u1=soln(i,1);
        u2=soln(i,2);
        u3=soln(i,3);
        u4=soln(i,4);
        energy_vec(i)=u2*u4-(u3)^2/2+(u2)^2+(1+params.mu)*u1^2/2-params.nu*u1^3/3+u1^4/4;
    end

    val=energy_vec;
   % Guess=E*norm(soln(1,:));
    %val=[Guess,energy_vec];
end