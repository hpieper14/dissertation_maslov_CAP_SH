function plots=plot_norm(phisol)
    N = max(size(phisol));
    norms=zeros(1,N);
    for i=1:N
        norms(i)=vecnorm(phisol(i,2:5));
    end
    figure
    plot(phisol(:,1),norms)
    xlabel('Time')
    ylabel('Euclidean Norm of Solution')
    title('Norm of Pulse Solution via Concatentation')
    plots=0;
end