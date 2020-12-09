function [ ModelParameter ] = Init_ModelParameter( model, env_nitrogen )

    %%% Initialize constraint vectors, and stoichiometric matrix
    vmin = model.lb;
    vmax = model.ub;
    Aeq = full(model.S);    % No external species exist in this model 
    beq = zeros(size(Aeq,1),1);

    for rxn = 1:length(vmin)
        if vmin(rxn) == -1000
            vmin(rxn) = -Inf;
        end
        if vmax(rxn) == 1000
            vmax(rxn) = Inf;
        end
    end
    

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Nutrients - general %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	vmin(strcmp(model.rxns,'Im_NH4')) = 0.0;    % --> NH4[c]
	vmax(strcmp(model.rxns,'Im_NH4')) = Inf;
	vmin(strcmp(model.rxns,'Im_NO3')) = 0.0;    % 2 ATP[c] + 2 H2O[c] --> NO3[c] + 2 ADP[c] + 2 Pi[c]
	vmax(strcmp(model.rxns,'Im_NO3')) = Inf;
	vmin(strcmp(model.rxns,'Im_Pi'))  = 0.0;    % 3 ATP[c] + 3 H2O[c] --> Pi[h] + 3 ADP[c] + 3 Pi[c]
	vmax(strcmp(model.rxns,'Im_Pi'))  = Inf;
	vmin(strcmp(model.rxns,'Im_SO4')) = 0.0;    % 3 ATP[c] + 3 H2O[c] --> SO4[c] + 3 ADP[c] + 3 Pi[c]
	vmax(strcmp(model.rxns,'Im_SO4')) = Inf;
	vmin(strcmp(model.rxns,'Im_H2S')) = 0.0;    % --> H2S[c]
	vmax(strcmp(model.rxns,'Im_H2S')) = Inf;
	vmin(strcmp(model.rxns,'Ex_O2'))  = -Inf;    % O2[c] -->
	vmax(strcmp(model.rxns,'Ex_O2'))  = Inf;
	vmin(strcmp(model.rxns,'Im_H2O')) = -Inf;    % --> H2O[c]
	vmax(strcmp(model.rxns,'Im_H2O')) = Inf;
	vmin(strmatch('Bio_',model.rxns)) = 0;	% Biomass functions
	vmax(strmatch('Bio_',model.rxns)) = 0;
	vmin(strcmp(model.rxns,'Im_CO2')) = -Inf;    % --> CO2[c]
	vmax(strcmp(model.rxns,'Im_CO2')) = Inf;

    
	%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Systemic changes %%%
	%%%%%%%%%%%%%%%%%%%%%%%%

	% vmin(strcmp(model.rxns,'Tr_G6P')) = 0.0;
	% vmax(strcmp(model.rxns,'Tr_G6P')) = Inf;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Nitrogen assimilation %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch env_nitrogen
		case 'both'
			a = 1; % dummy
        case 'nh4'
            vmax(strcmp(model.rxns,'Im_NO3')) = 0.0;
        case 'no3'
            vmax(strcmp(model.rxns,'Im_NH4')) = 0.0;
    end

    
    %%%%%%%%%%%%%
	%%% Other %%%
	%%%%%%%%%%%%%

    % Initialize exporter to carry zero flux despite oxygen efflux
    vmax(strmatch('Ex_', model.rxns)) = 0;
    vmin(strcmp(model.rxns,'Ex_O2'))  = -Inf;    % O2[c] -->
	vmax(strcmp(model.rxns,'Ex_O2'))  = Inf;

    % Initialize maintenance rxns to carry zero flux
    vmax(strmatch('NGAM_', model.rxns)) = 0;
    
    % Biomass reaction names
    biomRxnNames = {'Bio_AA'; 'Bio_CLim'; 'Bio_NLim'; 'Bio_opt'};
    
    ModelParameter.kNumbRxns = length(model.rxns);
	ModelParameter.kNumbMets = length(model.mets);
	
	ModelParameter.vmin = vmin;
	ModelParameter.vmax = vmax;
	ModelParameter.Aeq = Aeq;
	ModelParameter.beq = beq;
	ModelParameter.biomRxnNames = biomRxnNames;
    
end