function [weightsObjectiveDayspecific] = Read_weightTable()

    %%% read in weights table
    %
    %   - columns correspond to the timepoint(days) - 1st col is 0d here
    %   - 1st row corresponds to alpha
    %   - 2nd row corresponds to beta
    %   - 3rd row corresponds to gamma 
    
    % read in
    weightsObjectiveDayspecific = csvread('/Data/weightsObjectiveDayspecific.csv');

	
end