function [ ModelParameter ] = Impose_DiurnalModelParameters( ModelParameter, model, env_diurnal )
    
	vmin = ModelParameter.vmin;
	vmax = ModelParameter.vmax;
    
    %%% Energy source import IDs
	HNUtransID = find(strcmp(model.rxns,'Im_hnu'));	% --> hnu[h]
	GLCtransID = find(strcmp(model.rxns,'Ex_Glc'));	% Glc[c] -->
    
	if strcmp(env_diurnal, 'day') == 1
        env_trophic = 'auto';
		vmin(strcmp(model.rxns,'PPIF6PK_c')) = 0.0;   % Phosphofructokinase utilizing Pyrophosphate in compartment "c" (EC 2.7.1.90)
		vmax(strcmp(model.rxns,'PPIF6PK_c')) = 0.0;
		vmin(strcmp(model.rxns,'G6PDH_h')) = 0.0;   % Glucose-6-phosphate dehydrogenase in compartment "h" (EC 1.1.1.49)
		vmax(strcmp(model.rxns,'G6PDH_h')) = 0.0;
	elseif strcmp(env_diurnal, 'night') == 1
        env_trophic = 'hetero';
% 		vmin(strcmp(model.rxns,'RBC_h')) = 0.0;   % RuBisCO (utilizing CO2) in compartment "h" (EC 4.1.1.39)
% 		vmax(strcmp(model.rxns,'RBC_h')) = 0.0;
% 		vmin(strcmp(model.rxns,'RBO_h')) = 0.0;   % RuBisCO (utilizing O2) in compartment "h" (EC 4.1.1.39)
% 		vmax(strcmp(model.rxns,'RBO_h')) = 0.0;
% 
		vmin(strcmp(model.rxns,'GAPDH1_h')) = 0.0;   % Glyceraldehyde 3-phosphate dehydrogenase 1 in compartment "h" (EC 1.2.1.12)
		vmax(strcmp(model.rxns,'GAPDH1_h')) = 0.0;
		vmin(strcmp(model.rxns,'GAPDH2_h')) = 0.0;   % Glyceraldehyde 3-phosphate dehydrogenase 2 in compartment "h" (EC 1.2.1.12)
		vmax(strcmp(model.rxns,'GAPDH2_h')) = 0.0;

		vmin(strcmp(model.rxns,'FBPase_c')) = 0.0;   % Fructose-1,6-biphosphatase in compartment "c" (EC 3.1.3.11)
		vmax(strcmp(model.rxns,'FBPase_c')) = 0.0;
		vmin(strcmp(model.rxns,'FBPase_h')) = 0.0;   %  Fructose-1,6-biphosphatase in compartment "h"
		vmax(strcmp(model.rxns,'FBPase_h')) = 0.0;

		vmin(strcmp(model.rxns,'SBPase_h')) = 0.0;   % Sedoheptulosebiphosphatase in compartment "h" (EC 3.1.3.11 OR 3.1.3.37)
		vmax(strcmp(model.rxns,'SBPase_h')) = 0.0;

% 		vmin(strcmp(model.rxns,'Ru5PK_h')) = 0.0;   % Phosphoribulokinase in compartment "h" (EC 2.7.1.19)
% 		vmax(strcmp(model.rxns,'Ru5PK_h')) = 0.0;

		vmin(strcmp(model.rxns,'MalDH3_h')) = 0.0;   % Malate dehydrogenase (NADP+) in compartment "p" (EC 1.1.1.82)
		vmax(strcmp(model.rxns,'MalDH3_h')) = 0.0;

		vmin(strcmp(model.rxns,'ATPase_h')) = 0.0;
		vmax(strcmp(model.rxns,'ATPase_h')) = 0.0;
    end

    switch env_trophic
        case 'auto'
            % energy uptake
            vmin(HNUtransID) = 0.0; 
            vmax(HNUtransID) = 1000;
            vmin(GLCtransID) = 0.0;
            vmax(GLCtransID) = 0.0;
            % other
            vmin(strmatch('PSI_h',model.rxns)) = 0.0;
            vmax(strmatch('PSI_h',model.rxns)) = Inf;
            vmin(strmatch('PSII_h',model.rxns)) = 0.0;
            vmax(strmatch('PSII_h',model.rxns)) = Inf;
            vmin(strmatch('Fd-NADPR_h',model.rxns)) = 0.0;
            vmax(strmatch('Fd-NADPR_h',model.rxns)) = Inf;
            vmin(strcmp(model.rxns,'Im_CO2')) = -Inf;    % --> CO2[c]
            vmax(strcmp(model.rxns,'Im_CO2')) = Inf;
        case 'hetero'
            % energy uptake
            vmin(HNUtransID) = 0.0;
            vmax(HNUtransID) = 0.0;
            vmin(GLCtransID) = -1;
            vmax(GLCtransID) = 0.0;
            % other
            vmin(strmatch('PSI_h',model.rxns)) = 0.0;
            vmax(strmatch('PSI_h',model.rxns)) = 0.0;
            vmin(strmatch('PSII_h',model.rxns)) = 0.0;
            vmax(strmatch('PSII_h',model.rxns)) = 0.0;        
            vmin(strmatch('Fd-NADPR_h',model.rxns)) = -Inf;
            vmax(strmatch('Fd-NADPR_h',model.rxns)) = Inf;
            vmin(strcmp(model.rxns,'Im_CO2')) = -Inf;    % --> CO2[c]
            vmax(strcmp(model.rxns,'Im_CO2')) = 0.0;
    end
    
    ModelParameter.vmin = vmin;
	ModelParameter.vmax = vmax;
    
end