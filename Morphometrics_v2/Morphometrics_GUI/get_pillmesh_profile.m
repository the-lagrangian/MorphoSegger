%Tristan Ursell
%Pill Mesh Profiler
%July 2015
%

function sum1 = get_pillmesh_profile(X,Y,im_mat)

%get bounding box of polygon
x_l=floor(min(X));
x_r=ceil(max(X));

y_d=floor(min(Y));
y_u=ceil(max(Y));

sub_mat=im_mat(y_d:y_u,x_l:x_r);

%readjust coordinates
Xnew=X-(x_l-1);
Ynew=Y-(y_d-1);


%{
figure;
hold on
imagesc(sub_mat)
plot(X-(x_l-1),Y-(y_d-1),'k')
axis equal tight
%}

sum1=0;
for i=1:size(sub_mat,1)
    ypxl=[i-1/2,i+1/2,i+1/2,i-1/2];
    
    for j=1:size(sub_mat,2)
        xpxl=[j-1/2,j-1/2,j+1/2,j+1/2];
        
        
        %get pixel coords
        [x_int,y_int]=polybool('intersection',Xnew,Ynew,xpxl,ypxl); %Added the functions - Andres Florez 05/17/22
        %[x_int,y_int]=polybool_tsu(Xnew,Ynew,xpxl,ypxl);

        if ~isempty(x_int)
            %get intersection area
            area1=polyarea(x_int,y_int);
            
            sum1=sum1+sub_mat(i,j)*area1; 
            
            %plot(x_int,y_int,'k','linewidth',2)
            %pause(0.2)
        end
    end
end

%{
%test figure
figure;
hold on
imagesc(sub_mat)
plot(X-(x_l-1),Y-(y_d-1),'k')
axis equal tight

figure;
hold on
imagesc(im_mat)
plot(X,Y,'k')
axis equal tight
%}










