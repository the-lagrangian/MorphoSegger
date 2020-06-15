
%extract data from a cell array produced by analysis on morphometrics
function [foci,noFoci,cppd,dt]= extractFociData(dataName)

data_peaks_foci={};
data_no_foci_t={};
data_c_plus_d_t={};
data_doubling_t={};

fn = fieldnames(dataName);
for k=1:numel(fn)
    %if( isnumeric(dataName.(fn{k})) )
       data_peaks_foci{k,1} = dataName.(fn{k})(4,:);
       data_no_foci_t{k,1} = dataName.(fn{k})(6,:);
       data_c_plus_d_t{k,1} = dataName.(fn{k})(7,:);
       data_doubling_t{k,1} = dataName.(fn{k})(8,1);
    %end
end


foci=data_peaks_foci;
noFoci=data_no_foci_t;
cppd=data_c_plus_d_t;
dt=data_doubling_t;

%%
cell2csv("foci_data.csv",foci)

%save foci data to csv