function [ exp_name_matrix ] = Create_ExpNameVec( dataset_text )
%%CreateExpNameVec Create matrix with experiment details: 1st line -- name of mutant relative to -- 2nd line
    exp_name_matrix = dataset_text;
    exp_name_matrix = exp_name_matrix(1:3, 3:end);
    exp_names_ids   = (1:( (size(exp_name_matrix,2)+1) /4) ) * 4 - 3; % data element containing name is followed by three empty elements <- determine nonempty entries
    exp_name_matrix = exp_name_matrix([1, 3], exp_names_ids); % reduce to nonempty elements  
    
end
