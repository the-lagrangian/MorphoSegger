function [wavelet_bin] = stack_wavelet_foci(tsStack,noiseTh)

wavelet_bin=struct('data', repmat({zeros(size(tsStack(1).data))}, 1, size(tsStack,2)));
    
    for k=1:1:size(tsStack,2)
   
            A=tsStack(k).data;
            sigma=2;

            %A = histeq(A);
            A = imgaussfilt(A,sigma);

            %imshow(imadjust(B));

            bw=wavelet_foci2(A,noiseTh);

            %imshow(BW)



           %%
            % Compute the distance transform of the complement of the binary image.
            D1 = bwdist(~bw);
            %D = D(D>20);

            D1 = imgaussfilt(D1,2);
            % figure
            % imshow(D1,[],'InitialMagnification','fit')
            % title('Distance transform of ~bw')



            %%
            % Complement the distance transform, and force pixels that don't belong to
            % the objects to be at |-Inf| .
            %D=D(D<2);
            D1 = -D1;
            D1(~bw) = -Inf;
            %%


            % Compute the watershed transform and display the resulting label matrix as
            % an RGB image.
            L1 = watershed(D1);
        %    rgb1 = label2rgb(L1,'jet',[.5 .5 .5]);
        %    figure
        %    imshow(rgb1,'InitialMagnification','fit')
        %    title('Watershed transform of D')

            loc_B=find(L1<=1);

            L1(loc_B)=0;
            L1(~loc_B)=1;
            A(loc_B)=0;

            BW = imbinarize(L1);

            SE = strel('disk',1,0);
            BW = imerode(BW,SE);
            %SE = strel('disk',1,0);   
            BW = imdilate(BW,SE);
            
            wavelet_bin(k).data=BW;
            
    end      
end

