function [foci4,foci2,tb,tc,tau_cluster] = cell_cycle_main2(foci,foci3_N,kmeans_name)

foci_t=foci;
%tau_cluster=zeros(1,length(foci));

aux=find(foci_t>=0);  
                    
foci_t=foci_t(aux(1):aux(1)+5);
foci_t=foci_t(foci_t>0);
foci_t=mean(foci_t);
 
if foci_t<=1.5

[foci4,foci2,~,tb,tc] = cell_cycle2(foci,foci3_N);   
                                
else
    
[foci4,foci2,~,tb,tc] = cell_cycle(foci,foci3_N);
                               
end
                        
tau_cluster=clustering_kmeans(foci2);   
                   
savefig(kmeans_name); %save clustering figures

close all
                    
      
end                    
