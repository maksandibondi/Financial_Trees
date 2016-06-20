function [VaR, pdf, ddf] = VaR(path,num_of_iter,level)
%% VaR

% Profit/Loss variable is calculated as path(n,end)-path(1,1) that means
% difference between end and start price of the asset
step = 0.1;
pl_set = -path(1,1):step:2*path(1,1); % profit/loss variable scope
indic = 0; % indicator if VaR is already found
for i = 1:size(pl_set,2);
    counter_freq = 0;
    counter_density = 0;
    
    for k = 1:num_of_iter
        
        if i~=size(pl_set,2)
        
            if (path(k,end)-path(1,1))<pl_set(i+1) 
                counter_freq = counter_freq+1;
            end;
        
            if ((path(k,end)-path(1,1))<pl_set(i+1))&&((path(k,end)-path(1,1))>pl_set(i))
                counter_density = counter_density+1;
            end;
            
        else
            
            counter_freq = num_of_iter;
        
        end;
        
    end; 
    
    pdf(i) = counter_freq/num_of_iter;
    ddf(i) = counter_density/num_of_iter;
    
    if (pdf(i) >= 1-level)&&(indic==0)
        VaR = -pl_set(i);
        indic = indic + 1;
    end;
end;

return;

