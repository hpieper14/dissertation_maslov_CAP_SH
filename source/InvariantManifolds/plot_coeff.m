function plots=plot_coeff(coeff,order)
    normorder=zeros(order+1,order+1);
    plotpoints=[];
    suborder=1;
    while suborder<order+1
        for i=0:suborder
            for j = suborder-i
                vec=zeros(1,4);
                for k=1:4
                    vec(k)=coeff(i+1,j+1,k);
                end
                
                normpoint=norm(vec);
                normorder(i+1,j+1)=normpoint;
                point=[suborder;normpoint];
                plotpoints=[plotpoints,point];
            end
        end
        suborder=suborder+1;
    end
    
  plot(plotpoints(1,:), log(plotpoints(2,:)),'o');
  plots=0;
end