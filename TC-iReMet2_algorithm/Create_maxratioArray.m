function [ maxratio_array ] = Create_maxratioArray ( max_ratio)

    % init string array to save exponent numbers as strings in 
    maxratio_array = strings(1,max_ratio);

    % fill string array with exponents 
    for i=0:(max_ratio-1)   
        maxratio_array(i+1)= strcat( '1e',num2str(max_ratio-i));
    end
   
    % flip order for right order ( 54321 -> 12345)
    maxratio_array = fliplr(maxratio_array);
    % convert string to double so its useable as a number later on 
    maxratio_array = str2double(maxratio_array);
    
    
end 