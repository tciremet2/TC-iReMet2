function [ AddInformation ] = Init_informationarray()

% at day information the 
% 1col = all calculated maxratios
% 2col = all USED minratios
% 3col = all USED maxratios
%
%%% get rid of the hardcoded 10000 at some point 

% number of constraints made at timepoint
AddInformation.numberofconstraintsmade = ones(1,4)*NaN;
% table containing all calculated total maxratios for each timepoint
AddInformation.totalmaxratio = ones(10000,4)*NaN;
% table containg all used ratioconstraints min
AddInformation.usedminratio = ones(10000,4)*NaN;
% table containing all used ratioconstraints max
AddInformation.usedmaxratio = ones(10000,4)*NaN;

end

















