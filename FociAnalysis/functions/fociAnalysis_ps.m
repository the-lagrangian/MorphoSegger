function  [noiseTh] =  fociAnalysis_ps(stackname,frame,limits,paramFit,timeStep,Dparameter,exp_cut)

delta_D=5;

Dparameter=60;

N_D=1;

N_frames=20;

paramFit=N_frames;
%loading .tif stack and setting xy number
tsStack = tiffread(stackname);
xy_pos = char(regexp(stackname,'(_xy\w*).','match'));
xy_pos = [erase(xy_pos,"."),'_'];

th_levels=10;


posiciones=find(limits(:,2)==max(limits(:,2)) & limits(:,2)>max(limits(:,2))-N_frames+1);

if isempty(posiciones)==1
    
    warning("Error")
end    
u_frame=max(limits(:,2));
p_frame=u_frame-N_frames+1;

wavelet_bin_1=stack_wavelet_foci_ps(tsStack,th_levels,p_frame,u_frame,N_frames);

noiseTh=zeros(length(posiciones),1);

    
for z=1:length(posiciones)%Ncell
try     
     N=posiciones(z);

   
             
     l_cell=zeros(1,u_frame);
     foci=zeros(1,u_frame);
     foci2=zeros(1,u_frame); 
     foci3=zeros(1,u_frame); 
     foci4=zeros(1,u_frame); 
     tb=zeros(1,u_frame);
     tc=zeros(1,u_frame);
     a1=[];
     a2=[];
     a3=[];
     foci_m=zeros(N_D,u_frame);
     foci3_m=zeros(th_levels,u_frame);
     
     start=limits(N,2)-N_frames+1;
     finish=limits(N,2);
     cont=finish-start;
     finish=limits(N,2);

    

     %if cont<=paramFit-1  %consecutive points for a single cell to include the trajectory
     %     continue
     %else  
         %%% Build up the mother trajectory for kymograph. It uses 
         %%% the centerline of each cell for the last frame 

         %%%start


    
    %end
%         %%%ends
% 
%         %%%Generates and ROI using regionprops and BoundingBox
%         %%%in the last frame for each cell. This is
%         %%%used for wavelet foci counting. 
%         %%%start
% 
% 
% %    try
%         xs=round(frame(finish).object(ind).Xcont);
%         ys=round(frame(finish).object(ind).Ycont);
%         label=num2str(N);
% %    catch
% %        warning(label)
% %        continue
% %    end
% 
% 
%     AB=tsStack(start).data;
%     Im_s=zeros(size(AB));
% 
%     for p=1:length(xs)
%         Im_s(ys(p), xs(p))=1; 
%     end
% 
%     se = strel('square',3);
%     Im_s = imclose(Im_s,se);     
%     Im_s=imfill(Im_s,'holes');
%     stats = regionprops(Im_s,'BoundingBox');
% 
%     xMin = ceil(stats.BoundingBox(1));
%     xMax = xMin + stats.BoundingBox(3) - 1;
%     yMin = ceil(stats.BoundingBox(2));
%     yMax = yMin + stats.BoundingBox(4) - 1;
% 
%     % ends for for each frame each cell. 
%   
%       %Foci counting
        for k=start:finish      
        
        try    
            
            
            

               allCN = vertcat(frame(k).object.cellID); % comma separated list expansion 
               ind = find(allCN == N);
               

                xs=round(frame(k).object(ind).Xcont);
                ys=round(frame(k).object(ind).Ycont);
                
                
                l_cell(k)=frame(k).object(ind).length;   
                                
                A=tsStack(k).data;                                   
                AA=A;
                                                                                                                  
                A = imgaussfilt(A,2); 
                I=double(A);

                %%Generates binary per cell per frame
                Im1=zeros(size(I));

                for p=1:length(xs)
                    Im1(ys(p), xs(p))=1; 
                end

                se = strel('square',3);
                Im1 = imclose(Im1,se);     
                Im2=imfill(Im1,'holes');
                loc=find(Im2==0);
                loc1=find(Im2==1);

                % Generate the 3 trayectories based on the centerline (and
                % 2 parallel trajectories). The signals is integrated for
                % thress trajectories. It is referenced to the main
                % trajectory (last frame) using centroid.

                x=round(frame(k).object(ind).centerline(:,2));
                y=round(frame(k).object(ind).centerline(:,1));
                x_m=round(frame(k).object(ind).pill_mesh(:,2));
                y_m=round(frame(k).object(ind).pill_mesh(:,1));

                %x_ms2=round(frame(k).object(ind).pill_mesh(:,4));
                %y_ms2=round(frame(k).object(ind).pill_mesh(:,3));

                delta_x= round(abs(x-x_m)/2.5);
                delta_y= round(abs(y-y_m)/2.5);
                %[x,y]=cleanup(x,y);
                %l_cell(k)=length(x);

                lin_ind = sub2ind(size(I),x,y);
                lin_ind2 = sub2ind(size(I),x+delta_x,y+delta_y);
                lin_ind3 = sub2ind(size(I),x-delta_x,y-delta_y);

                a1=double(I(lin_ind));
                a2=double(I(lin_ind2));
                a3=double(I(lin_ind3));
                a=a1+a2+a3;
                aa=A(lin_ind);
                
                %some cells are problematic, use this to prevent stopping
                %the calculation
%                 try
                b=sgolayfilt(a,4,11); 
%                 catch
%                    continue
%                 end

                [pks1,locs1]=findpeaks(b,'MinPeakDistance',5);
                

                % DIEGO'S FOCI COUNTING: Since singal2noise is low, local
                % variability of each peak was assessed by the derivate of the
                % signal over 6 neighboors pixels. Lower than Diego_parameter is
                % rejected.
                
                for p=1:N_D
                                       
                contar=0;
                
                for q=1:length(pks1)
                      try
                           bbb=(b(locs1(q)-3:locs1(q)+3));
                      catch
                      continue
                      end
                      
                      condition_f=sum(abs(diff(bbb)));
  
                      if condition_f>(Dparameter+(p-1)*delta_D)
                         contar=contar+1;    
                      end
                end

                foci_m(p,k)=contar;
               
                end
      
               %%%% Wavelet counting
               
               for w=1:th_levels
               
               
               AB=wavelet_bin_1(w,k).data;
      
               AB(loc)=0;
                                   
               stats = regionprops(AB,'Area');
              
               stats2=vertcat(stats.Area);
                                   
               stats2=stats2(stats2>15);
                                   
               foci3_m(w,k)=length(stats2);
                                   
               end
               
               %%%% End wavelet 
             
               
        
        catch        
            continue
        end       
               
               
        end

    aux2=l_cell==0;
    l_cell(aux2)=NaN;         
    
    exp_cut=10;
    paramFit=10;
    time=timeStep*(0:1:length(l_cell)-1); 
    [~, ~, t_t]=createFit2_exp(time, l_cell,paramFit,exp_cut);

    if t_t==1      
        continue                                
    end    
    
    close all

    
    [t_n]=cell_cycle_opt_ps(foci_m,foci3_m);
    
    noiseTh(z)=t_n;
    

catch        
    continue
end     
end



end

    