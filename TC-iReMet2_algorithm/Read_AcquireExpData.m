function [ ratio_min, ratio_max, exp_description ] = Read_AcquireExpData ( norm_name, kNumbSDratioConstr, experimentname )

	% load data
    eval( strcat( 'data_read = Read_mixed_csv(''/Data/MTWT_ratio_mets_',experimentname,'_mean_sdev_sorted_fitforreadin_', norm_name, '.csv'', ''\t'', 70); ' ) );
   
    % remove the columns without values
	% odd columns: mean relative changes
	% even columns: relative standard errors
	met_ratios_values = Create_NumericMatrix( data_read );
	
	% Experiment details: 1st row -- mutant, relative to -- 2nd row
	exp_description = Create_ExpNameVec( data_read );

	% Determine min- and max-ratio w.r.t. error intervals

	error_ids   = [1 : ( size(met_ratios_values, 2) / 2 )] * 2;
	mean_ids = error_ids - 1;
	ratio_mean = met_ratios_values( :, mean_ids );
	ratio_error = met_ratios_values( :, error_ids );
	ratio_min = ratio_mean - kNumbSDratioConstr * ratio_error;
	ratio_max = ratio_mean + kNumbSDratioConstr * ratio_error;

    
    
    ratio_min(ratio_min < 0) = 0;  % concentration ratio cannot be smaller than zero
	
end