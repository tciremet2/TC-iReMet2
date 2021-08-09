function [ measuredModelMetIDs_dataIDs, met_const_ratio_ids ] = Read_GeneralData ( model, experimentname )

	% Load (presumably) constant metabolite names -- a problem with a string in the file occurs, however, the needed data is correctly imported -- TODO: fix
    [ met_const_ratio, dummy1, dummy2 ] = xlsread('/Data/mets_constant_arnold_model.xls', 1);   
    met_const_ratio = met_const_ratio(:, 1);
	met_const_ratio_ids = find(met_const_ratio == 1);  % indices of the metabolites which are assumed to be constant

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Create mapping from measured metabolite indices found in the model to metabolite indices in the dataset
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	% read mapping of measured metabolite full name indices in the dataset and
	% metabolite full name indices in the model
    eval( strcat( '[ dummy1, metNames_expVSmodel, dummy2 ] = xlsread(''/Data/MTWT_ratio_mets_dictionary_',experimentname,'_sorted_model_data_dictionary_corrected.xls'',1);' ) );
    
	% create 1. cell array with measured metabolite names found in the model,
	% 2., vector with corresponding indices in the dataset
	%%% 3:end to 2:end cause header only one line in own case 
    metNames_expVSmodel = metNames_expVSmodel(2:end, :); % delete 1st rows containing headers; column 1: metabolite names (dataset), column 2: metabolite names (model)
	measuredMetNames_inModel = {};
    measuredMetIDs_inExp = [];
	for row = 1 : size(metNames_expVSmodel, 1)
		if ~isequal(metNames_expVSmodel{row, 2}, '')
			measuredMetNames_inModel = [measuredMetNames_inModel; metNames_expVSmodel{row, 2}];
			measuredMetIDs_inExp = [measuredMetIDs_inExp; row];
		end
	end

	% create matrix with metabolite indices in the model corresponding to found
	% metabolite names and the corresponding metabolite indices in the dataset
	measuredModelMetIDs_dataIDs = []; % column 1: ids of measured metabolites in the model, column 2: ids of measured metabolites in the dataset
	for metNameID = 1 : length(model.metNames)
		search_id = find(strcmp(measuredMetNames_inModel, model.metNames(metNameID)));
		if isempty(search_id) ~= 1
			measuredModelMetIDs_dataIDs = [ measuredModelMetIDs_dataIDs; 0, 0 ];
			measuredModelMetIDs_dataIDs(end, 1) = metNameID;
			measuredModelMetIDs_dataIDs(end, 2) = measuredMetIDs_inExp(search_id);
		end
    end
	
end