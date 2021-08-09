function [ fopt, hnuOpt ] = Determine_BiomassProdAndPhotonInflux( model, vmin, vmax, biomass_reaction, knockout )

    biomass_reaction_id       = find(strcmp(model.rxns, biomass_reaction));
    vmax(biomass_reaction_id) = Inf;    % Set upper biomass boundary

    % tomlab FBA
    c_bm                      = zeros(length(model.rxns), 1);
    c_bm(biomass_reaction_id) = 1;
    v                         = tom('v',length(model.rxns), 1);
    
    % use this part to set fluxes of knockoutmutants to 0
    %{
    if ~isempty(strfind(knockout, 'glyk'))                 
        vmax(strcmp(model.rxns, 'PGAK_h')) = 0.0;
        vmin(strcmp(model.rxns, 'PGAK_h')) = 0.0;
    end
    if ~isempty(strfind(knockout, 'shm')) 
        vmax(strcmp(model.rxns, 'GlyHMT_m')) = 0.0;
        vmin(strcmp(model.rxns, 'GlyHMT_m')) = 0.0;
    end
    if ~isempty(strfind(knockout, 'hpr')) 
        vmax(strcmp(model.rxns, 'GCEADH_p')) = 0.0;
        vmin(strcmp(model.rxns, 'GCEADH_p')) = 0.0;
    end
    if ~isempty(strfind(knockout, 'pglp')) 
        vmax(strcmp(model.rxns, 'PGP_h')) = 0.0;
        vmin(strcmp(model.rxns, 'PGP_h')) = 0.0;
    end
    %}
    
    objective   = -1*c_bm'*v;
    constraints = {model.S*v==0; vmin<=v<=vmax};
    Prob        = sym2prob(objective, constraints); % Prob=sym2prob('',objective,constraints);
    Result      = tomRun('snopt', Prob,1);
    fopt        = -1*Result.f_k;  % maximum biomass production
    hnuOpt      = Result.x_k(strcmp('Im_hnu', model.rxns));  % corresponding photon import

end
