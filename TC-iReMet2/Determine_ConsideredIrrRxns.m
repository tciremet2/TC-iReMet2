function [ considered_irrev_ids_vec, artificial_rxns_and_transporters ] = Determine_ConsideredIrrRxns( model, vmin )

    % Identify considered reactions: irrversible, no transporters
    all_irrev_ids_vec                = find(vmin == 0);
    artificial_rxns_and_transporters = [strmatch('NGAM_', model.rxns),  % maintenance reactions
                                strmatch('Bio_', model.rxns),  % biomass reactions
                                strmatch('Im_', model.rxns),  % import reactions
                                strmatch('Ex_', model.rxns),  % export reactions
                                strmatch('Tr_', model.rxns),  % compartment transport reactions
                                find(strcmp(model.rxns, 'Si_H'))   % proton source/sink
                                ];
    considered_irrev_ids_vec = setdiff(1:length(model.rxns), artificial_rxns_and_transporters);
    considered_irrev_ids_vec = intersect(considered_irrev_ids_vec, all_irrev_ids_vec);
    
end