function [] = Evaluate_results (maxratio,toppercent, kNumbExps, model ) 

    % read in data 
    eval(strcat('data=dlmread(''/home/mpimp-golm.mpg.de/pries1347/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/day/without_constCofactorRatios/with_slacks/with_slack_minimization/1e',num2str(maxratio),'/rep1/deviationMatrix_-_norm_to_Col0.dat',''');'));
    
    % get knumbrxns
    kNumbRxns = size(model.S,2);
    % init vars
    ToppercentRxnNumber = ceil(kNumbRxns*toppercent);
    
    % calc top % for every timepoint, transport rxns excluded
    for timepoint = 1:kNumbExps
        eval(strcat('toprxnsTime',num2str(timepoint),'= Determine_topRxns(data(:,',num2str(timepoint),'), ToppercentRxnNumber, model);'));
        %disp(strcat('toprxnsTime',num2str(timepoint),'= Determine_topRxns(data(:,',num2str(timepoint),'), ToppercentRxnNumber, model);'));
    end
    
    % here conserved thingy 
    %toprxnsTime5 = [];    
    
    % top rxns indecies
    topbig = [];
    for i=1:kNumbExps+1
        eval(strcat('topbig = [ topbig, toprxnsTime',num2str(i),'];'));
    end
        
    % new stuff here 
    a = intersect(topbig(:,1),intersect(topbig(:,4),intersect(topbig(:,7),topbig(:,10))));
    b = intersect(topbig(:,2),intersect(topbig(:,5),intersect(topbig(:,8),topbig(:,11))));
    c = intersect(topbig(:,3),intersect(topbig(:,6),intersect(topbig(:,9),topbig(:,12))));
    
    % fill the table with the respective rxn names 
    big = cell(size(topbig,1),(kNumbExps+1)*3);
    for i=1:size(topbig,1)
        for j=1:size(topbig,2)
            big(i,j)= model.rxns(topbig(i,j));
        end
    end
    
    % header names
    header = {'Rxns_WTMT_flux_0d','Rxns_MTWT_flux_0d','overallRxns_absOfWTMT_flux_0d',...
                'Rxns_WTMT_flux_1d','Rxns_MTWT_flux_1d','overallRxns_absOfWTMT_flux_1d',...
                'Rxns_WTMT_flux_7d','Rxns_MTWT_flux_7d','overallRxns_absOfWTMT_flux_7d',...
                'Rxns_WTMT_flux_21d','Rxns_MTWT_flux_21d','overallRxns_absOfWTMT_flux_21d',...
                'dummy1','dummy2','dummy3'};
    
    % add empty cols for comfort 
    topbig = [topbig, zeros(size(topbig,1),3)];          
            
    writetable(cell2table(big,'VariableNames',header),'/home/mpimp-golm.mpg.de/pries1347/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/toprxns.csv','Delimiter','\t');
    writetable(array2table(topbig,'VariableNames',header),'/home/mpimp-golm.mpg.de/pries1347/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/toprxnIDs.csv','Delimiter','\t');
            
                

end 