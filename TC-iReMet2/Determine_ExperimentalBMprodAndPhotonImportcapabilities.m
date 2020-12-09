function [ fopt, hnuOpt ] = Determine_ExperimentalBMprodAndPhotonImportcapabilities( model, ModelParameter, knockoutMutants )    

    fopt   = NaN * ones(length(ModelParameter.biomRxnNames), length(knockoutMutants));
    hnuOpt = NaN * ones(length(ModelParameter.biomRxnNames), length(knockoutMutants));

    for biomRxn = 1:length(ModelParameter.biomRxnNames)
        for koM = 1:length(knockoutMutants)
            vmin = ModelParameter.vmin;
            vmax = ModelParameter.vmax;
            %if koM == 2  % special case shm mutant     % use this if you need to adjust specific fluxes in order for the mutants to produce biomass
            %    vmax(483:487) = Inf;
            %end
	
            [ biomass_production, photon_import ] = Determine_BiomassProdAndPhotonInflux( model, vmin, vmax, ModelParameter.biomRxnNames{biomRxn}, knockoutMutants{koM} );
            
            if biomass_production > 1e-5
                fopt(biomRxn, koM)   = biomass_production; 
                hnuOpt(biomRxn, koM) = photon_import;
            else
                fopt(biomRxn, koM)   = 0;
                hnuOpt(biomRxn, koM) = 0;
            end
        end
    end
    
end