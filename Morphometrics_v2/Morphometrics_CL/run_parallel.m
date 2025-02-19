% runMorphometrics                              %
% Input: Analysis folder, parameter file 
%
% This function runs morphometrics_mask_cl on multiple xy points in parallel
% Andres Florez - 04/29/2020
% andrewflorez@gmail.com
% Harvard University

function run_parallel(dirname, params)

%dirname=pwd;
%dirname=fixDir(dirname);
cleanMorphometrics(dirname); % delete temporary files before run
contents = dir([dirname,'xy*']); %List all xy folders
num_dir_tmp = numel(contents);
nxy = [];
num_xy = 0;

%Creates a struct for all seg directories 
dirnamelist=cell(1,num_dir_tmp);
for i = 1:num_dir_tmp
    if (contents(i).isdir) && (numel(contents(i).name) > 2)
        num_xy = num_xy+1;
        nxy = [nxy, str2double(contents(i).name(3:end))];
        dirnamelist{i} = [dirname,contents(i).name,filesep,'seg'];
    end
end

%List seg mask files and runs morphometrics in parallel

parfor (k = 1:num_xy)
    dirname_xy = dirnamelist{k}; %added
    mask = dir([dirname_xy,filesep,'mask1seg_xy*.tif']);
    maskFile_tmp= [mask(1).folder,filesep, mask(1).name]; %important to select the mask file and prevent calling Gsegt etc.
    morphometrics_mask_cl(maskFile_tmp,params,[],1); 
end

% move morphometrics files to independent folder
for q = 1:num_xy
    movefile([dirnamelist{q},filesep,'*_CONTOURS*.mat'], [dirname,contents(q).name,filesep,'morphometrics']);
end

disp('Parallel Morphometrics done.')
clearvars -except dirname
% %workers=6;
% if  workers % shutting down parallel pool     
%     delete(poolobj);
% end

end