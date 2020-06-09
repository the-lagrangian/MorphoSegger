function [total,dist] = dist_clustering_kmeans(foci)


foci_n=foci(foci>0);

x=[1:1:length(foci_n)];

aux=foci_n>=2 & foci_n<=10;

foci_n=foci_n(aux);

x=x(aux);

pos_a=diff(x)<=2;

x=x(pos_a);
    
foci_n=foci_n(pos_a);

if length(foci_n)<10
    
    return
    
end    

dist=0;

total=0;


XXX=[x',foci_n'];

try
[idx,~] = kmeans(XXX,3,'MaxIter',100,'Replicates',10);%,'Start','sample');%,'Distance','cosine');'Distance','cosine'


for i=1:3
    
    a_vec=XXX(idx==i,2);
   
    total=total+length(a_vec);
    dist=dist+sum(diff(a_vec).^2);%/length(a_vec);
    
    
end
catch
return    
end

end

