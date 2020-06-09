function [foci,foci3] = cell_cycle_opt(foci_m,foci3_m)

opt_matrix=zeros(size(foci3_m,1),size(foci_m,1));


for i=1:size(foci3_m,1)
    
    
    for j=1:size(foci_m,1)


            foci_t=foci_m(j,:);

            aux=find(foci_t>=0);  

            foci_t=foci_t(aux(1):aux(1)+5);
            foci_t=foci_t(foci_t>0);
            foci_t=mean(foci_t);

                                if foci_t<=1.5

                                            [~,~,~,tb,tc] = cell_cycle2(foci_m(j,:),(foci_m(j,:)+foci3_m(i,:))/2);   



                                else

                                            [~,~,~,tb,tc] = cell_cycle(foci_m(j,:),(foci_m(j,:)+foci3_m(i,:))/2);


                                end
                                
                                [total, dist_cl]=dist_clustering_kmeans(foci_m(j,:));                                 
                      
                                dist_b=((tb(3)-tb(2))^2+(tb(4)-tb(3))^2);
                             
                                dist_c=((tc(3)-tc(2))^2+(tc(4)-tc(3))^2);
                                
                                opt_matrix(i,j)=(tb(4)+tb(3)+tb(2))+(tc(3)+tc(2))+total-(dist_cl+dist_b+dist_c);
                                
                                %opt_matrix(i,j)=dist_cl;

    end                   
                            
end            
 
loc_max=find(opt_matrix==max(max(opt_matrix)));

[row,col] = ind2sub(size(opt_matrix),loc_max(1));

foci3=foci3_m(row,:);

foci=foci_m(col,:);

      
end                    
