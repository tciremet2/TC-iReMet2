function [ OptResults ] = Init_QFVA( AnalysisParameters, OptResults, artificial_rxns_and_transporters, model )
    
    AP = AnalysisParameters;
    
    kNumbRxns = size(OptResults.deviationMatrix, 1);
    kNumbExps = size(OptResults.deviationMatrix, 2);
    
    considered_rxns = 1:kNumbRxns;
    if AP.qpFVA_with_artificial_rxns == 0
        considered_rxns = setdiff(considered_rxns, artificial_rxns_and_transporters);
    end
    
    considered_experiments  = intersect( find(~isnan(OptResults.deviations)), find(sqrt(OptResults.deviations) > 0.001) );
    mean_absolut_deviations = mean( abs( OptResults.deviationMatrix( considered_rxns, considered_experiments ) ), 2);
    [ ~, sorted_ids ]       = sort(mean_absolut_deviations, 1, 'descend');
    number_of_qpFVA_rxns    = floor( AP.top_rxn_ratio_qpFVA * length(considered_rxns) );
    examined_rxns           = considered_rxns( sorted_ids(1:number_of_qpFVA_rxns) );
    
    % Determine size of parts
    average_exps_per_part          = ceil( length(considered_experiments) / AP.qpFVA_kNumbParts );
    qpFVA_parts_of_considered_exps = [];
    for part = 1:(AP.qpFVA_kNumbParts-1)
        qpFVA_parts_of_considered_exps = [qpFVA_parts_of_considered_exps, average_exps_per_part];
    end
    qpFVA_parts_of_considered_exps = [qpFVA_parts_of_considered_exps, length(considered_experiments) - (AP.qpFVA_kNumbParts-1)*average_exps_per_part];
    
    OptResults.qpFVA_examined_rxns            = examined_rxns;
    OptResults.qpFVA_examined_rxnIDs          = model.rxns( examined_rxns );
    OptResults.qpFVA_parts_of_considered_exps = qpFVA_parts_of_considered_exps;
    OptResults.qpFVA_considered_experiments   = considered_experiments;
    
    OptResults.fluxRangeMatrixWT       = ones(length(examined_rxns), 2*kNumbExps)*NaN;
    OptResults.fluxRangeMatrixM        = ones(length(examined_rxns), 2*kNumbExps)*NaN;
    OptResults.fluxRangeDiffMatrixWT   = ones(length(examined_rxns), kNumbExps)*NaN;
    OptResults.fluxRangeDiffMatrixM    = ones(length(examined_rxns), kNumbExps)*NaN;
    
end
