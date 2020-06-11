
function MakeResultsMatrix()
% This function deletes Gparent and Gsegt (temporary) files when
% Morphometrics stops due to errors.

% Andres Florez - 04/29/2020
% andrewflorez@gmail.com
% Harvard University


% Select folder with mat files cc_results

disp('Select Folder with ND2 files');
dirname = uigetdir();
dirname = fixDir(dirname);








end