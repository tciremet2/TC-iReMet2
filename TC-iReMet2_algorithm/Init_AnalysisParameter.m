function [ AnalysisParameters ] = Init_AnalysisParameter( filename )

    fid = fopen( filename, 'r' );
    tline = fgetl(fid);
    while ischar(tline)
        eval(tline)
        tline = fgetl(fid);
    end
    fclose(fid);

    AnalysisParameters.env_nitrogen               = env_nitrogen;
    AnalysisParameters.kNumbSDratioConstr         = kNumbSDratioConstr;  % ratio constraint boundaries: mean +/- kNumbSDratioConstr * standard error
    AnalysisParameters.min_mutant_BM_prod         = min_mutant_BM_prod;  % minimum biomass production of the mutant relative to the wildtype
    AnalysisParameters.with_slacks                = with_slacks;  % options: '1' (yes) and '0' (no)
    AnalysisParameters.with_slack_minimization    = with_slack_minimization;  % options: '1' (yes) and '0' (no)
    AnalysisParameters.slack_weighting            = slack_weighting;  % with low weighting, optmization may lead to an increased flux distance but decreased sum of absolute ratio constraint slacks
    AnalysisParameters.with_constCofactorRatios   = with_constCofactorRatios;  % options: '1' (yes) and '0' (no)
    AnalysisParameters.kNumbTestX_0               = kNumbTestX_0;  % number of trials of different random initial solutions
    AnalysisParameters.max_ratio                  = max_ratio;  % largest total ratio constraint; limited to prevent badly scaled problem
    AnalysisParameters.allowed_rel_dev_linear     = allowed_rel_dev_linear;  % allowed relative deviation of individual linear constraints in qFVA
    AnalysisParameters.allowed_rel_dev_nonlinear  = allowed_rel_dev_nonlinear;  % allowed relative deviation of individual nonlinear constraints in qFVA
    %%% 
    AnalysisParameters.experimentname             = experimentname; % name of the experiment to run code on 
    AnalysisParameters.with_rangeBiom             = with_rangeBiom;
    AnalysisParameters.devrangeBiom               = devrangeBiom;
    
    if AnalysisParameters.with_slack_minimization == 1
        AnalysisParameters.kNumbSlacksPerRatioConstr = 1;
    else
        AnalysisParameters.kNumbSlacksPerRatioConstr = 1;
    end
    
    AnalysisParameters.kNumbRepetition            = kNumbRepetition; % number of reruns for each cutoff-value for evaluation of the stability of the results  
    AnalysisParameters.with_daySpecificObjective  = with_daySpecificObjective;
    
    
end

