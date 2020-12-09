function [ ConstraintInfo ] = Determine_ConstraintInfos( model, measuredModelMetIDs_dataIDs, met_const_ratio_ids, considered_irrev_ids_vec, Tconstrained_rxns )
  
    % create index vector of considered reactions with at least one measured substrate ratio
	reduced_stoich_mat = model.S( measuredModelMetIDs_dataIDs(:, 1), :);
	constrained_rxns = find( any( reduced_stoich_mat < 0, 1 ) ); % indices of reactions with at least one measured substrate ratio
    % find idx of reactions that can be constrained by either Metsratio or Tratio
    constrained_rxns = union(Tconstrained_rxns,constrained_rxns); % indecies of reactions with either measured Metratio or Tratio
	constrained_rxns = intersect( constrained_rxns, considered_irrev_ids_vec ); % reduce to rxns that are among the considered rxns

	% create vector of considered metabolites 
	reduced_stoich_mat = model.S( :, constrained_rxns );
	ratio_constraint_mets = find( any(reduced_stoich_mat < 0, 2 ) );

	% create index vectors for ratio_constraint_mets identifying measured and unmeasured
	% metabolites, as well as indicators of their model ids
	constraint_mets_measured   = find( ismember( ratio_constraint_mets, measuredModelMetIDs_dataIDs(:, 1) ) );
	constraint_mets_unmeasured = find( ~ismember( ratio_constraint_mets, measuredModelMetIDs_dataIDs(:, 1) ) );
	constraint_mets_unmeasured_constant = intersect( constraint_mets_unmeasured, met_const_ratio_ids);

    ConstraintInfo.constrained_rxns = constrained_rxns;
	ConstraintInfo.kNumbRatioConstr = length(constrained_rxns);
	ConstraintInfo.constraint_mets = ratio_constraint_mets;
	ConstraintInfo.constraint_mets_measured = constraint_mets_measured;
	ConstraintInfo.constraint_mets_unmeasured = constraint_mets_unmeasured;
	ConstraintInfo.constraint_mets_unmeasured_constant = constraint_mets_unmeasured_constant;
	ConstraintInfo.measuredModelMetIDs_dataIDs = measuredModelMetIDs_dataIDs;
	ConstraintInfo.kNumbRatioConstr = length(ConstraintInfo.constrained_rxns);
    ConstraintInfo.Tconstrained_rxns = Tconstrained_rxns;

end
