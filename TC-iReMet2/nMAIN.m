function [] = nMAIN()    
    
    %%%%%%%%%%%%%%%%%%
    %%%  run vars  %%%
    %%%%%%%%%%%%%%%%%%
    
    InitFilename = 'AnalysisParameters.init';
    env_diurnal = 'day';
    %cutoff = 8;
    %repetition =1;
    %day = 1;

    
    %%%%%%%%%%%%%%%%%% 
    %%% Load model %%% 
    %%%%%%%%%%%%%%%%%% 
    
    load('Model_A');
    model = Model_A; 
    model = Correct_model_mistakes(model);
    
    
    %%%%%%%%%%%%%%%%%%
    %%% Initialize %%%
    %%%%%%%%%%%%%%%%%%
    
    AnalysisParameters  = Init_AnalysisParameter( InitFilename );    
    resultPath          = Determine_ResultPath( AnalysisParameters, env_diurnal );
    ModelParameterInit  = Init_ModelParameter( model, AnalysisParameters.env_nitrogen );
    

    %%%%%%%%%%%%%%%%%
    %%% Read Data %%%
    %%%%%%%%%%%%%%%%%
    
    norm_names = {'norm_to_Col0'};
    for norm_id = 1:length(norm_names)
        [ Mratio_min, Mratio_max, experiment_description ] = Read_AcquireExpData( norm_names{norm_id}, AnalysisParameters.kNumbSDratioConstr, AnalysisParameters.experimentname );
        ExpData.(norm_names{norm_id}).Mratio_min = Mratio_min;
        ExpData.(norm_names{norm_id}).Mratio_max = Mratio_max;
        ExpData.(norm_names{norm_id}).description = experiment_description;
    end   
    [ Tratio_min, Tratio_max, Tconstrained_rxns] = Read_transcriptomicData( AnalysisParameters.kNumbSDratioConstr, AnalysisParameters.experimentname); % transcriptomic data
    [ measuredModelMetIDs_dataIDs, met_const_ratio_ids ] = Read_GeneralData( model, AnalysisParameters.experimentname );
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Determine model parameter and constraint information depending on day/night %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ModelParameter = Impose_DiurnalModelParameters( ModelParameterInit, model, env_diurnal);
    [ constrRxnsIrr, artificial_rxns_and_transporters ] = Determine_ConsideredIrrRxns( model, ModelParameter.vmin );
    ConstraintInfo = Determine_ConstraintInfos( model, measuredModelMetIDs_dataIDs, met_const_ratio_ids, constrRxnsIrr, Tconstrained_rxns );  
    fileID = fopen(strcat(resultPath, 'Ratio_constrained_irrev_rxns.dat'),'w');
    for col = 1:length(ConstraintInfo.constrained_rxns)
        fprintf( fileID, '%s\t%s\n', model.rxns{ConstraintInfo.constrained_rxns(col)},  model.rxnNames{ConstraintInfo.constrained_rxns(col)} );
    end
    fclose(fileID);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Effects of KO on biomass and photon import %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    KOeffects.knockoutMutants = {'noKO'}; % for knockout use reactionname to knockout, e.g. {'glyk', 'shm', 'hpr', 'pglp', 'noKO'}
    [ KOeffects.fopt, KOeffects.hnuOpt ] = Determine_ExperimentalBMprodAndPhotonImportcapabilities( model, ModelParameter, KOeffects.knockoutMutants );

    
    %%%%%%%%%%%%%%%%%%%%
    %%% Optimization %%%
    %%%%%%%%%%%%%%%%%%%%
    
    OptInit = Init_OptArrays( model, ModelParameter.kNumbRxns, ModelParameter.kNumbMets ); 
    
    for norm_id = 1:length(norm_names)
        
        % set general parameters
        kNumbRxns = size(model.S, 2);
        considered_experiments = Determine_ConsideredExperiments( ExpData.(norm_names{norm_id}), env_diurnal, AnalysisParameters.experimentname );
        kNumbExps = length(considered_experiments);
        OptResults = Init_OptResults(kNumbExps, kNumbRxns, ConstraintInfo.kNumbRatioConstr );
        maxratio_array = Create_maxratioArray( AnalysisParameters.max_ratio);
        rvCorrStruct = Init_rvCorrStruct(maxratio_array,AnalysisParameters,kNumbRxns,kNumbExps);
        weightsObjectiveDayspecific = Read_weightTable();
        
        for cutoff = 1:length(maxratio_array)
            
            % set maxratio
            AnalysisParameters.max_ratio = maxratio_array(cutoff);
            
            for repetition = 1:AnalysisParameters.kNumbRepetition
               
                % set resultpath to corresponding maxratio directory                
                resultPathDirectory = Set_resultPathToDirectory(resultPath,cutoff,repetition);  
                fluxArrayDayBefore  = Init_fluxArrayDayBefore(kNumbRxns);
                % init information array
                AddInformation = Init_informationarray();
                
                for day = 1:length(considered_experiments)
                    
                                      
                    % show information for better understanding in command line
                    disp(considered_experiments(day))
                    disp(ExpData.(norm_names{norm_id}).description{1, considered_experiments(day)})
                    fprintf('cutoff: %d \n',cutoff)
                    fprintf('repetition: %d \n',repetition)
                    
                    % set metabolomic and transcriptomic min-/max ratio 
                    % according to their dayspecific values
                    measured_met_ratio_min = ExpData.(norm_names{norm_id}).Mratio_min(:, considered_experiments(day));
                    measured_met_ratio_max = ExpData.(norm_names{norm_id}).Mratio_max(:, considered_experiments(day));
                    measured_T_ratio_min   = Tratio_min(:,day);
                    measured_T_ratio_max   = Tratio_max(:,day);
                    
                    [ ExpCondBoundaries ] = Impose_ExperimentSpecificBoundaries( ModelParameter, KOeffects, ExpData.(norm_names{norm_id}), considered_experiments(day), model, AnalysisParameters );
                    
                    [RatioConstraint,AddInformation] = Define_RxnRatioConstraintsWslacks( model, measured_met_ratio_min, measured_met_ratio_max, measured_T_ratio_min, measured_T_ratio_max, ConstraintInfo, AnalysisParameters, day, AddInformation);

                    OptInput = Assemble_OptInput (AnalysisParameters, OptInit, ConstraintInfo.kNumbRatioConstr, ...
                                                  ExpCondBoundaries, RatioConstraint, kNumbRxns, fluxArrayDayBefore, day, weightsObjectiveDayspecific );                                                                                                              

                    Result = Run_TomlabWithDifferentStartingValues ( OptInit.Name, OptInput, AnalysisParameters.kNumbTestX_0 );
                    
                    ResultsPlusConstraints.(strcat('experiment', num2str(day))).Result = Result;    % Store results and constraints for qpFVA
                    ResultsPlusConstraints.(strcat('experiment', num2str(day))).OptInput = OptInput;
                    ResultsPlusConstraints.(strcat('experiment', num2str(day))).ExpCondBoundaries = ExpCondBoundaries;
                    
                    OptResults = Store_OptResults( OptResults, Result, OptInput, day, kNumbRxns, AnalysisParameters, ConstraintInfo, RatioConstraint );
                    % array set for vminvmax von distances
                    fluxArrayDayBefore = Set_fluxArrayDayBefore(fluxArrayDayBefore, OptResults, day);
                                              
                end
                
                Write_ResultsToFiles( OptResults, AnalysisParameters, ExpData, norm_names{norm_id}, InitFilename, resultPathDirectory, considered_experiments, ResultsPlusConstraints, kNumbExps );
                Write_Additionalinformation(AddInformation,AnalysisParameters,resultPathDirectory,cutoff);
                rvCorrStruct = Store_WTandMT_fluxMatrix ( rvCorrStruct, OptResults, cutoff, repetition, kNumbExps); 
                
                if AnalysisParameters.with_qpFVA == 1
                    qFVAstoreFileName = strcat(resultPathDirectory, 'QFVArequiredData_', norm_names{norm_id}, '.mat');
                    OptResults        = Init_QFVA( AnalysisParameters, OptResults, artificial_rxns_and_transporters, model );
                    save( qFVAstoreFileName, 'ResultsPlusConstraints', 'AnalysisParameters', 'OptResults' );
                end
                
                % no qFVA since sampling was done
                
            end          
            
        end
        
        % plot correlation of cutoffs 
        Plot_correlationOfCutoffs(rvCorrStruct,AnalysisParameters,kNumbExps,maxratio_array)
        
        % evaluate pathways tables according to different measures -- don't
        % use server for this!!! => only run local
        %{
        % do heatmaps of pathways, relative and nominal changes, etc. 
        Eval_MAIN();
        % do enrichment analysis of pathways
        Eval_enrichment();
        %}
        
    end
end
