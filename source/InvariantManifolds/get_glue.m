function [bestL, index] = get_glue(Lsol, Rsol, time_vec)
    diff_vec=Lsol-Rsol;
   
    figure
    
    if min(size(diff_vec))==1
        plot(time_vec, abs(diff_vec))
        [val, index] = min(abs(diff_vec));
    else
        plot(time_vec, vecnorm(diff_vec.'))
        [val, index] = min(vecnorm(diff_vec.'));
    end
    
    
    
    
    bestL=time_vec(index);
    disp('Optimal time of gluing:')
    disp(bestL);
    disp('Size of gap at optimal glue location:')
    disp(val);
    
  
    
end
