function [ resultPath ] = Determine_ResultPath( AnalysisParameters, env_diurnal )

    if AnalysisParameters.with_slacks == 1
        if AnalysisParameters.with_slack_minimization == 1
            slackUtilizationPath = 'with_slacks/with_slack_minimization/';
        else
            slackUtilizationPath = 'with_slacks/without_slack_minimization/';
        end
    else
        slackUtilizationPath = 'without_slacks/';
    end

    if AnalysisParameters.with_constCofactorRatios == 1
        resultPath = strcat('Results/', env_diurnal, '/with_constCofactorRatios/');
    else
        resultPath = strcat('Results/', env_diurnal, '/without_constCofactorRatios/');
    end

    resultPath = strcat(resultPath, slackUtilizationPath);

end
