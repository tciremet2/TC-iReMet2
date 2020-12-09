function [] = Write_ResultsToFiles ( OptResults, AnalysisParameters, ExpData, experimentIDCond, InitFilename, resultPath, considered_experiments, ResultsPlusConstraints, kNumbExps )

    fileID = fopen( strcat(resultPath, 'Deviations_-_', experimentIDCond, '.dat'), 'w' );
    for experimentID = 1:length(OptResults.deviations)
        experimentIDName = [ExpData.(experimentIDCond).description{1, considered_experiments(experimentID)}, ' - ', ExpData.(experimentIDCond).description{2, considered_experiments(experimentID)}];
        fprintf(fileID, '%s\t', experimentIDName);
        fprintf(fileID, '%s\n', sqrt(OptResults.deviations(experimentID)));
    end
    fclose(fileID);

    dlmwrite(strcat(resultPath, 'deviationMatrix_-_', experimentIDCond, '.dat'), OptResults.deviationMatrix, '\t');
    dlmwrite(strcat(resultPath, 'wildtypeFluxMatrix_-_', experimentIDCond, '.dat'), OptResults.wildtypeFluxMatrix, '\t');
    dlmwrite(strcat(resultPath, 'mutantFluxMatrix_-_', experimentIDCond, '.dat'), OptResults.mutantFluxMatrix, '\t');

    dlmwrite(strcat(resultPath, 'fluxRangeMatrixWT_-_', experimentIDCond, '.dat'), OptResults.fluxRangeMatrixWT, '\t');
    dlmwrite(strcat(resultPath, 'fluxRangeMatrixM_-_', experimentIDCond, '.dat'), OptResults.fluxRangeMatrixM, '\t');

    dlmwrite(strcat(resultPath, 'fluxRangeDiffMatrixWT_-_', experimentIDCond, '.dat'), OptResults.fluxRangeDiffMatrixWT, '\t');
    dlmwrite(strcat(resultPath, 'fluxRangeDiffMatrixM_-_', experimentIDCond, '.dat'), OptResults.fluxRangeDiffMatrixM, '\t');

    dlmwrite(strcat(resultPath, 'slackValueMatrix_-_', experimentIDCond, '.dat'), OptResults.slackValueMatrix, '\t');

    dlmwrite(strcat(resultPath, 'ratioConstraintMatrix_-_', experimentIDCond, '.dat'), OptResults.ratioConstraintMatrix, '\t');

    copyfile( InitFilename, strcat(resultPath, 'AnalysisParameters.init') );
    
    if AnalysisParameters.with_qpFVA == 1
        fileID = fopen( strcat(resultPath, 'Rxns_with_largest_deviation', experimentIDCond, '.dat'), 'w' );
        for rxn = 1:length(OptResults.qpFVA_examined_rxnIDs)
            fprintf(fileID, '%s\n', OptResults.qpFVA_examined_rxnIDs{rxn});
        end
        fclose(fileID);
    end

    %{
    % save the linear part of the objective in table for qfva as it cant be
    % given as a parameter to the dc function in lpconassign
    linConPart = ones(length(ResultsPlusConstraints.experiment1.OptInput.d),kNumbExps)*NaN;
    for exp = 1:kNumbExps
       eval(strcat( 'linConPart(:,',num2str(exp),') = ResultsPlusConstraints.experiment',num2str(exp),'.OptInput.d;') );
    end
    dlmwrite(strcat(resultPath, 'linear_objectivePart_experiments.dat'),linConPart,'\t'); 
    
    % save const addition part of the objective to extra file 
    constadd = ones(1,kNumbExps)*NaN;
    for exp = 1:kNumbExps        
        eval(strcat('constadd(' ,num2str(exp), ')=ResultsPlusConstraints.experiment' ,num2str(exp), '.OptInput.const_add;'));
    end
    dlmwrite(strcat(resultPath, 'const_add_objective.dat'),constadd,'\t')
    
    % test for x_k
    X0 = ones(length(ResultsPlusConstraints.experiment1.Result.x_k),kNumbExps)*NaN;
    for exp = 1:kNumbExps        
        eval(strcat('X0(:,' ,num2str(exp), ')=ResultsPlusConstraints.experiment' ,num2str(exp), '.Result.x_k;'));
    end
    dlmwrite(strcat(resultPath, 'Xnull.dat'),X0,'\t')
    %}

end