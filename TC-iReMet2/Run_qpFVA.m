function [ ] = Run_qpFVA( InitFilename, env_diurnal, partID )

    AP         = Init_AnalysisParameter( InitFilename );
    resultPath = Set_resultPathToDirectory(resultPath,cutoff,repetition);

    disp(strcat('You''re performing qpFVA part ', blanks(2), num2str(partID), '/', num2str(AP.qpFVA_kNumbParts)));
    
    load(strcat(resultPath, 'QFVArequiredData_norm_to_Col0.mat'));
    
    kNumbPrecedingExps = sum(OptResults.qpFVA_parts_of_considered_exps(1:(partID-1)));
        
    for experiment = OptResults.qpFVA_considered_experiments( (kNumbPrecedingExps+1) : (kNumbPrecedingExps + OptResults.qpFVA_parts_of_considered_exps(partID)) )'
    %for experiment = 1:4
        
        Result   = ResultsPlusConstraints.(strcat('experiment', num2str(experiment))).Result;
        OptInput = ResultsPlusConstraints.(strcat('experiment', num2str(experiment))).OptInput;
        ExpCondBoundaries = ResultsPlusConstraints.(strcat('experiment', num2str(experiment))).ExpCondBoundaries;
        
        [ fluxRangeWT, fluxRangeM ] = Run_IndividualQFVA( Result, OptResults, ExpCondBoundaries, OptInput, AP, experiment );

        OptResults.fluxRangeMatrixWT(:, [(2*experiment-1), 2*experiment]) = fluxRangeWT;
        OptResults.fluxRangeMatrixM(:, [(2*experiment-1), 2*experiment])  = fluxRangeM;
        OptResults.fluxRangeDiffMatrixWT(:, experiment) = fluxRangeWT(:, 2) - fluxRangeWT(:, 1);
        OptResults.fluxRangeDiffMatrixM(:, experiment)  = fluxRangeM(:, 2) - fluxRangeM (:, 1);
        
    end  
        
    OptResults_currentID = strcat('OptResults_qFVA_part_', num2str(partID), '_of_', num2str(AP.qpFVA_kNumbParts));
    eval(strcat(OptResults_currentID, ' = OptResults'));
    save( strcat(resultPath, OptResults_currentID, '.mat'), OptResults_currentID );

end