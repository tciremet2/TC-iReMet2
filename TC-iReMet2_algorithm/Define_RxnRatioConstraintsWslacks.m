function  [RatioConstraint,AddInformation] = Define_RxnRatioConstraintsWslacks( model, measured_met_ratio_min, measured_met_ratio_max, measured_T_ratio_min, measured_T_ratio_max,...
                                                                                ConstraintInfo, AnalysisParameters,day, AddInformation)
% constraint_matrix has as many rows as linear constraints derived from
% measured metabolite abundances and twice as many columns as reactions
    
    CI = ConstraintInfo;
    AP = AnalysisParameters;
    kNumbRxns = size(model.S, 2);
    % Initialize constraint cell array
    
    if AP.with_slacks == 0
        constraint_matrix = zeros(2*length(CI.constrained_rxns), (2*kNumbRxns));
    else
        constraint_matrix = zeros(2*length(CI.constrained_rxns), (2*kNumbRxns + AP.kNumbSlacksPerRatioConstr * length(CI.constrained_rxns)));
    end
    
    constraint_matrix_for_readout = zeros(length(CI.constrained_rxns), 2);
    
    % Initialize min/max ratios
    allExpM_ratio_min = min(measured_met_ratio_min);
    allExpM_ratio_max = max(measured_met_ratio_max);
    allExpT_ratio_min = min(measured_T_ratio_min);  
    allExpT_ratio_max = max(measured_T_ratio_max);  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Reaction flux constraints %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    id_in_mapping_matrix = find( ismember (CI.measuredModelMetIDs_dataIDs(:, 1), CI.constraint_mets(CI.constraint_mets_measured) ) );
    
    ratio_min = zeros(length(CI.constraint_mets), 1);
    ratio_max = zeros(length(CI.constraint_mets), 1);
    
    
	ratio_min( CI.constraint_mets_measured ) = measured_met_ratio_min( CI.measuredModelMetIDs_dataIDs(id_in_mapping_matrix , 2) );
	ratio_max( CI.constraint_mets_measured ) = measured_met_ratio_max( CI.measuredModelMetIDs_dataIDs(id_in_mapping_matrix , 2) );
	ratio_min( CI.constraint_mets_unmeasured ) = ones( length( CI.constraint_mets_unmeasured ), 1) * allExpM_ratio_min; % CAUTION: here using minimum of all measured ratios - I may also use the minimum of ratios used in the optimization
	ratio_max( CI.constraint_mets_unmeasured ) = ones( length( CI.constraint_mets_unmeasured ), 1) * allExpM_ratio_max;
    if AP.with_constCofactorRatios == 1
        ratio_min( CI.constraint_mets_unmeasured_constant ) = ones(length(CI.constraint_mets_unmeasured_constant), 1); 
        ratio_max( CI.constraint_mets_unmeasured_constant ) = ones(length(CI.constraint_mets_unmeasured_constant), 1);
    end
    
    %%%  information variables
    numberOfMadeConstraints = 0;
    totalmaxratiotable = [];
    usedmaxratiotable = [];
    usedminratiotable = []; 
    for id = 1 : length(CI.constrained_rxns) % id = 15
        rxn = CI.constrained_rxns(id);
        % Determine substrates of rxn
        substrate_ids = find(model.S(:, rxn) < 0);
        % the absolute of the substrates' stochiometric coefficients is
        % required

        total_ratio_max = prod( power( ratio_max( ismember(CI.constraint_mets, substrate_ids) ), abs( model.S(substrate_ids, rxn) ) ) );       
%  here add transcript ratio once
        if isequaln(measured_T_ratio_max(rxn),NaN)
            total_ratio_max = total_ratio_max*(1*allExpT_ratio_max);  %  if Tratio is unmeasured 
        else
            total_ratio_max = total_ratio_max*(1*measured_T_ratio_max(rxn));  %  if Tratio is measured
        end
        %total_ratio_max = total_ratio_max*(1*measured_T_ratio_max(rxn));  %  HOW TO SET DEVIATION MULTIPLIER
        totalmaxratiotable = [totalmaxratiotable; total_ratio_max];
        
        if total_ratio_max < AP.max_ratio  % To prevent numerical instabilites, constraints with too large values are neglected
            constraint_matrix(2*id - 1, rxn)           = -1 * total_ratio_max;
            constraint_matrix(2*id - 1, rxn+kNumbRxns) = 1;
            constraint_matrix_for_readout(id, 2)       = total_ratio_max; 
            if AP.with_slacks == 1
                % constraint_matrix(2*id - 1, 2*kNumbRxns + id) = 1;  %
                % slack entry, upper part original, changed to -1 so the
                % slack value gets used properly 
                constraint_matrix(2*id - 1, 2*kNumbRxns + id) = -1;  % slack entry
            end
        
            total_ratio_min = prod( power( ratio_min( ismember(CI.constraint_mets, substrate_ids) ), abs( model.S(substrate_ids, rxn) ) ) );
            %  
            if isequaln(measured_T_ratio_min(rxn),NaN)
                total_ratio_min = total_ratio_min*(1*allExpT_ratio_min);  % if Tratio is unmeasured 
            else
                total_ratio_min = total_ratio_min*(1*measured_T_ratio_min(rxn));  % if Tratio is measured
            end
            %total_ratio_min = total_ratio_min*(1*measured_T_ratio_min(rxn)); %  

            constraint_matrix(2*id, rxn)           = -1 * total_ratio_min;
            constraint_matrix(2*id, rxn+kNumbRxns) = 1;
            constraint_matrix_for_readout(id, 1)   = total_ratio_min;
            if AP.with_slacks == 1
                constraint_matrix(2*id , 2*kNumbRxns + id)  = 1;    % slack entry
            end 
            
            %%% information variable 
            numberOfMadeConstraints = numberOfMadeConstraints+1;
            usedmaxratiotable = [usedmaxratiotable; total_ratio_max];
            usedminratiotable = [usedminratiotable; total_ratio_min];
        end
    end

    
    ratio_b_L = zeros( 2 * ConstraintInfo.kNumbRatioConstr, 1 );
    ratio_b_L( (1:ConstraintInfo.kNumbRatioConstr) * 2 - 1) = -Inf;

    ratio_b_U = zeros( 2 * ConstraintInfo.kNumbRatioConstr, 1 );
    ratio_b_U( (1:ConstraintInfo.kNumbRatioConstr) * 2) = Inf;
    
    RatioConstraint.ratio_constraint_matrix       = constraint_matrix;
    RatioConstraint.constraint_matrix_for_readout = constraint_matrix_for_readout;
    RatioConstraint.ratio_b_L                     = ratio_b_L;
    RatioConstraint.ratio_b_U                     = ratio_b_U;
    
    %%% information object 
    AddInformation.numberofconstraintsmade(1,day) = numberOfMadeConstraints;
    AddInformation.totalmaxratio(1:length(totalmaxratiotable),day) = totalmaxratiotable;
    AddInformation.usedminratio(1:length(usedminratiotable),day) = usedminratiotable;
    AddInformation.usedmaxratio(1:length(usedmaxratiotable),day) = usedmaxratiotable;
    
    %writetable(totalmaxratiotable,'_totalmaxratiotable'); %  
end