function [] = My_qfva (InitFilename, env_diurnal, resultPath)

    % loading needed parts 
    AP = Init_AnalysisParameter(InitFilename);
    load(strcat(resultPath, 'QFVArequiredData_norm_to_Col0.mat'));
    %kNumbExps = size(OptResults.deviations,1);
    kNumbExps = 4;
    
    % run QFVA for all experiments
    for experiment = 1:kNumbExps
       
        % get experimentspecific parts
        Result   = ResultsPlusConstraints.(strcat('experiment', num2str(experiment))).Result;
        OptInput = ResultsPlusConstraints.(strcat('experiment', num2str(experiment))).OptInput;
        ExpCondBoundaries = ResultsPlusConstraints.(strcat('experiment', num2str(experiment))).ExpCondBoundaries;
        save('optintputest.mat','OptInput');
        
        % run qfa 
        [ fluxRangeWT, fluxRangeM ] = Run_IndividualQFVA( Result, OptResults, ExpCondBoundaries, OptInput, AP, experiment );
        
        % save results of qfva
        OptResults.fluxRangeMatrixWT(:, [(2*experiment-1), 2*experiment]) = fluxRangeWT;
        OptResults.fluxRangeMatrixM(:, [(2*experiment-1), 2*experiment])  = fluxRangeM;
        OptResults.fluxRangeDiffMatrixWT(:, experiment) = fluxRangeWT(:, 2) - fluxRangeWT(:, 1);
        OptResults.fluxRangeDiffMatrixM(:, experiment)  = fluxRangeM(:, 2) - fluxRangeM(:, 1);
                
    end
    
    % save results like this for now 
    save('test.mat');

end