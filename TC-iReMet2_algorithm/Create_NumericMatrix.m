function [ numericData ] = Create_NumericMatrix( input )

    numericData = input(:, 3:end); %# cut out first two columns
    numericData = numericData(4:end, :); %# cut out first three rows

    numbericDataIDs = 2 * [1:(size(numericData, 2)+1)/2] - 1;   %# ids of odd column entries
    numericData = numericData(:, numbericDataIDs);              %# cut out even column entries
    numericData = str2double(numericData);                      %# convert to double

end
