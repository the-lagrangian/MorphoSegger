function MakeResultsMatrix()
%Generates clist files from filaments data to be processed with
%gateTool from SuperSegger

% Andres Florez - 04/29/2020
% andrewflorez@gmail.com
% Harvard University


% Select folder with mat files cc_results

disp('Select Folder cc_res files');
dirname = uigetdir();
dirname = fixDir(dirname);

% Get the .mat files
contents=dir([dirname '*.mat']);
num_im = numel(contents);

%Calculate growth rate:
%growth_rate = (log(len_1) - log(len_0)) ./ age;

data = [];
def={};


for i = 1:num_im
     data_c = loaderInternal([dirname,contents(i).name]);
     
     data(i,1) = data_c.output(2,1); %doubling time
%      data(i,2) = data_c.output(7,2); %B period 2 foci
%      data(i,3) = data_c.output(7,3); %B period 4 foci
%      data(i,4) = data_c.output(7,4); %B period 8 foci
%      
%      %Diego's method
%      data(i,5) = data_c.output(8,2); %C period 2 foci %Dudoso
%      data(i,7) = data_c.output(8,3); %C period 4 foci
%      data(i,9) = data_c.output(8,4); %C period 8 foci
%      
%      %clustering
%      data(i,6) = data_c.output(9,2); %C period 2 foci 
%      data(i,8) = data_c.output(9,3); %C period 4 foci
%      data(i,10) = data_c.output(9,4); %C period 8 foci
     
     [setter]  = clistSetter ();
     def = setter(:)';
     
     
     %rows: time-dependent 2D array
     %1. length
     %3. Pole 1 variation
     %4. Pole 2 variation
     %5. Foci time course wavelet
     %6. Foci time Diego's method
     save('clist','data','def')

end
end

function [setter] = clistSetter ()
%Variable names

setter = [{'Doubling time'};
          {'B period 2nd foci'};
          {'B period 4nd foci'};
          {'B period 8nd foci'};
          {'C period 2nd foci (warning)'};
          {'C period 4nd foci'};
          {'C period 8nd foci'};
          {'C period 2nd foci (clustering)'};
          {'C period 4nd foci (clustering)'};
          {'C period 8nd foci (clustering)'};
          ];

end

function data = loaderInternal(filename)
data = load( filename );
end
