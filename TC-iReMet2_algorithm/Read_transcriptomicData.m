function [Tratio_min, Tratio_max, Tconstrained_rxns] = Read_transcriptomicData(kNumbSDratioConstr, experimentname)

    %%% load data  
    eval( strcat( 'tdata = readtable(''/Data/MTWT_ratio_transcriptomicdata_', experimentname,'_mean_dev.txt'',''Delimiter'',''\t'');' ) );
    
    %%% seperate ratio and deviation % ugly hard code, change this
    Tratio = table2array(tdata(:,[2 4 6 8])); % transcriptomic ratio 
    Tdev   = table2array(tdata(:,[3 5 7 9])); % transcriptomic ratio deviation  
    
    % true for measured reactions with measured transcription ratio 
    Tconstrained_rxns = find(all(~isnan(Tratio),2));
        
    % Determine min- and max-ratio w.r.t. error intervals
    Tratio_min = Tratio - kNumbSDratioConstr * Tdev;
	Tratio_max = Tratio + kNumbSDratioConstr * Tdev;
    
    Tratio_min(Tratio_min < 0) = 0;  % Transcriptomics ratio cant be smaller than 0
 
end