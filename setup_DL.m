%Add Paths for DL experiments

%NOTE: some of the comparison codes require MEX files to be built. Compiled
%mex code for 64-bit MATLAB running on MAC OS X are included. These files
%will need to be built for other platforms


%Get absolute path to the folder containing this file
basePath = [fileparts(mfilename('fullpath')) filesep];

addpath([basePath '/BiGAMP']) %BiG-AMP code
addpath([basePath '/TruboCS']) %TruboCS code



