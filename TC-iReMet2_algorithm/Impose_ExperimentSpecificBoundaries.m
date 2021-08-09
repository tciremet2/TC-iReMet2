function [ ExpCondBoundaries ] = Impose_ExperimentSpecificBoundaries( ModelParameter, KOeffects, ExpDataSpecific, day, model, AnalysisParameters )

    kNumbRxns = length(model.rxns);
        
    % Experiment-specific biomass reactions
    % MUTANT
    %%% changed everything so the outcome always is 'Bio_opt'
    if ~isempty(strfind(ExpDataSpecific.description{1, day}, AnalysisParameters.experimentname)) 
		biomRxnM   = 'Bio_opt';
		biomRxnIdM = find(strcmp(model.rxns, 'Bio_opt'));
	else
		biomRxnM   = 'Bio_opt'; % changed just to bioopt to be save
		biomRxnIdM = find(strcmp(model.rxns, 'Bio_opt')); % changed to bioopt to be save 
    end
    % WILDTYPE
    if ~isempty(strfind(ExpDataSpecific.description{2, day}, AnalysisParameters.experimentname))
		biomRxnWT   = 'Bio_opt';
		biomRxnIdWT = find(strcmp(model.rxns, 'Bio_opt'));
	else
		biomRxnWT   = 'Bio_opt'; % changed just to bioopt to be save
		biomRxnIdWT = find(strcmp(model.rxns, 'Bio_opt')); % changed to bioopt to be save 
    end             
	
%	knockoutRxns = {'PGAK_h'; 'GlyHMT_m'; 'GCEADH_p'; 'PGP_h'};  % order matches order of KOeffects.knockoutMutants
	
    % Initialize flux boundaries
    vminWT = ModelParameter.vmin; % WildType
    vmaxWT = ModelParameter.vmax;
    vminM = ModelParameter.vmin; % Mutant
    vmaxM = ModelParameter.vmax;
	
    foptM  = KOeffects.fopt( strcmp(ModelParameter.biomRxnNames, biomRxnM), strcmp(KOeffects.knockoutMutants, 'noKO') );  % KOeffects.knockoutMutants = {'glyk', 'shm', 'hpr', 'pglp', 'noKO'}, ModelParameter.biomRxnNames = {'Bio_AA'; 'Bio_CLim'; 'Bio_NLim'; 'Bio_opt'}
    foptWT = KOeffects.fopt( strcmp(ModelParameter.biomRxnNames, biomRxnWT), strcmp(KOeffects.knockoutMutants, 'noKO') );
    
%     vmaxM(strcmp(model.rxns, 'Im_hnu')) = KOeffects.hnuOpt( strcmp(ModelParameter.biomRxnNames, biomRxnM), strcmp(KOeffects.knockoutMutants, 'noKO') );
%     vmaxWT(strcmp(model.rxns, 'Im_hnu')) = KOeffects.hnuOpt( strcmp(ModelParameter.biomRxnNames, biomRxnWT), strcmp(KOeffects.knockoutMutants, 'noKO') );
%{
    for koID = 1:length(KOeffects.knockoutMutants)
    %%% MUTANT SPECIFIC BOUNDARIES
    if ~isempty(strfind(ExpDataSpecific.description{1, exp}, KOeffects.knockoutMutants{koID}))
        vminM(strcmp(model.rxns, knockoutRxns{koID})) = 0.0;
        vmaxM(strcmp(model.rxns, knockoutRxns{koID})) = 0.0;
        vmaxM(strcmp(model.rxns, 'Im_hnu')) = KOeffects.hnuOpt( strcmp(ModelParameter.biomRxnNames, biomRxnM), strcmp(KOeffects.knockoutMutants, KOeffects.knockoutMutants{koID}) );
        foptM = KOeffects.fopt( strcmp(ModelParameter.biomRxnNames, biomRxnM), strcmp(KOeffects.knockoutMutants, KOeffects.knockoutMutants{koID}) ); 
    end
    %%% WILDTYPE SPECIFIC BOUNDARIES
    if ~isempty(strfind(ExpDataSpecific.description{2, exp}, KOeffects.knockoutMutants{koID}))
        vminWT(strcmp(model.rxns, knockoutRxns{koID})) = 0.0;
        vmaxWT(strcmp(model.rxns, knockoutRxns{koID})) = 0.0;
        vmaxWT(strcmp(model.rxns, 'Im_hnu')) = KOeffects.hnuOpt( strcmp(ModelParameter.biomRxnNames, biomRxnWT), strcmp(KOeffects.knockoutMutants, KOeffects.knockoutMutants{koID}) );
        foptWT = KOeffects.fopt( strcmp(ModelParameter.biomRxnNames, biomRxnWT), strcmp(KOeffects.knockoutMutants, KOeffects.knockoutMutants{koID}) ); 
    end
    end

	if ~isempty(strfind(ExpDataSpecific.description{1, exp}, 'shm'))
		vmaxM([483:487])=Inf; % Allow amino acid efflux since otherwise shm mutant cannot produce biomass
	end
	if ~isempty(strfind(ExpDataSpecific.description{2, exp}, 'shm')) 
		vmaxWT([483:487])=Inf; % Allow amino acid efflux since otherwise shm mutant cannot produce biomass
    end
%}
    
    %%% Set biom constraints and percentage of biomproduction according to
    %%% specified approach in Analysisparameters.init
    if isequal(AnalysisParameters.with_rangeBiom,1)
       
        % read in biomass fraction table with time point specific biomass fraction values
        biomassFraction_time = table2array(readtable('/Data/percentageOfBiomassProductionAtEachTimepoint2.csv'));

        dummymutantminBiom = biomassFraction_time(1,day);
        if AnalysisParameters.devrangeBiom > 0 
            vminM(biomRxnIdM) = (dummymutantminBiom-AnalysisParameters.devrangeBiom)*foptM;
            if (dummymutantminBiom+AnalysisParameters.devrangeBiom) > 1
                vmaxM(biomRxnIdM) = foptM;
            else 
                vmaxM(biomRxnIdM) = (dummymutantminBiom+AnalysisParameters.devrangeBiom)*foptM;
            end
            vminWT(biomRxnIdWT) = foptWT;
            vmaxWT(biomRxnIdWT) = foptWT;
        else
            vminWT(biomRxnIdWT) = foptWT;
            vmaxWT(biomRxnIdWT) = foptWT;
            vminM(biomRxnIdM) = dummymutantminBiom * foptM;
            vmaxM(biomRxnIdM) = foptM; 
        end
    else
        
        vminWT(biomRxnIdWT) = foptWT;
        vmaxWT(biomRxnIdWT) = foptWT;
        vminM(biomRxnIdM) = AnalysisParameters.min_mutant_BM_prod * foptM;
        vmaxM(biomRxnIdM) = foptM;
        
    end
   
    ExpCondBoundaries.vminWT = vminWT;
    ExpCondBoundaries.vmaxWT = vmaxWT; 
    ExpCondBoundaries.vminM = vminM;
    ExpCondBoundaries.vmaxM = vmaxM;
    ExpCondBoundaries.biomRxnIdWT = biomRxnIdWT;
    ExpCondBoundaries.biomRxnIdM = biomRxnIdM;
        
end