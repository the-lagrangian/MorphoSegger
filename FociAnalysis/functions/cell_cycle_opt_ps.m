function [noiseTh] = cell_cycle_opt_ps(foci_m,foci3_m)

opt_matrix=zeros(size(foci3_m,1),size(foci_m,1));


for i=1:size(foci3_m,1)
    
    
    for j=1:size(foci_m,1)

    opt_matrix(i,j)=sum(foci3_m(j,:))^2-sum((foci_m(j,:)-foci3_m(i,:)).^2);
            
    end                   
                            
end            
 
loc_max=find(opt_matrix==max(max(opt_matrix)));

[row,~] = ind2sub(size(opt_matrix),loc_max(1));

noiseTh=row+4;


      
end                    
