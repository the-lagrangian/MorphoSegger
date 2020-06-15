
function [clist]=MakeResultsMatrix()



% Andres Florez - 04/29/2020
% andrewflorez@gmail.com
% Harvard University

% First Row: cell length.
% 2: Fitting: Doubling time and Rsqr. 
% 3: differences frame by frame of pole up
% 4: differences frame by frame of pole down\
% 5: Diego's foci counting 
% 6: Wavelet foci counting 
% 7: B times:1foci,2foci,4foci,8foci,16foci,32foci
% 8: C Times:1foci,2foci,4foci,8foci,16foci,32foci
% 9: Synchonization of foci counitng. Check
%                    first coumn
%                     -1: starts Ctime:4 
%                     -2: starts Btime:2 
%                     -3: starts Ctime:2   
%                     -4: starts Btime:1 

% Select folder with mat files cc_results

disp('Select Folder cc_res files');
dirname = uigetdir();
dirname = fixDir(dirname);

% Get the .mat files
contents=dir([dirname '*.mat']);
num_im = numel(contents);

%Calculate growth rate:
%growth_rate = (log(len_1) - log(len_0)) ./ age;

clist.data = [];
clist.def={};

for i = 1:num_im
     data_c = loaderInternal([dirname,contents(i).name]);
     clist.data(i,1) = data_c.output(2,1);
end



end



function data = loaderInternal(filename)
data = load( filename );
end
