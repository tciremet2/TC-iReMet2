function [OV] = Plot_findAndPlotMostRegulatedPathwaysAndReactionsOverTimepoints()


    %%%%%%%%%%%%
    %%% INIT %%%
    %%%%%%%%%%%%

    %%% read in tables
    % read in rxnlist with corresponding pathways 
    rxnpathway = readtable('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Data/RxnsList_ID_shortname_pathway.csv');
    % read in deviationmatrix of fluxes between WTMT
    deviationMatrix = dlmread('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/day/without_constCofactorRatios/with_slacks/with_slack_minimization/1e8/rep1/deviationMatrix_-_norm_to_Col0.dat','\t');
    % norm to max each col
    for col=1:size(deviationMatrix,2)
        deviationMatrix(:,col) = deviationMatrix(:,col)/max(abs(deviationMatrix(:,col)));
    end
    
    %%% 
    % convert to cell for easier handling 
    rxnpathway = table2cell(rxnpathway);
    % extract pathways and exclude empty cells 
    pathways = rxnpathway(:,6);
    pathways = pathways(~cellfun('isempty',pathways));
    
    %%%
    % create array where 1 corresponds to present reaction in pathway(column) 
    foundpathway = ones(size(rxnpathway,1),size(pathways,1))*NaN;
    for pathway = 1:size(pathways,1)
        foundpathway(:,pathway)=contains(rxnpathway(:,4),pathways{pathway});
    end
       
    %%%
    % create table of fluxdeviation and reactionIDs for plotting 
    for pathway = 1:size(foundpathway,2)
        eval(strcat ('pway.p' ,num2str(pathway), ' = [deviationMatrix(find(foundpathway(:,' ,num2str(pathway), ')),:), find(foundpathway(:,' ,num2str(pathway), '))];'));
    end
    
    
    %%% Init struc containing all important information and arrays
    OV.rxnpathway = rxnpathway;
    OV.deviationMatrix = deviationMatrix;
    OV.pathways = pathways; 
    OV.foundpathway = foundpathway;
    OV.overalltimepoints = [];
    OV.t12 = [];
    OV.t23 = [];
    OV.t34 = [];
    
    
    
    
    %%%%%%%%%%%%%%%%%%%
    %%% CALCULATION %%%
    %%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%% OVER ALL TIMEPOINTS    
    %%% pathways
    % vector of overall regulation of each pathway in order to rank them
    regulationPathways = ones(size(pathways,1),1)*NaN;
    % calculate overall regulation of pathway     
    for pathway = 1:size(pathways,1)
        % assign pathwayspecific distancearray
        eval(strcat('x = pway.p',num2str(pathway),';'));
        % calculate distances for timepoints
        distance1 = sum(abs(x(:,1) - x(:,2)));
        distance2 = sum(abs(x(:,2) - x(:,3)));
        distance3 = sum(abs(x(:,3) - x(:,4)));
        overalldist = distance1 + distance2 + distance3;
        % assign pathwayspecific overall regulation value
        regulationPathways(pathway) = overalldist;
    end
    % rank them 
    [sortvaluesPathway,pathwayID] = sort(regulationPathways,'descend');
    regulatedPathwaysSorted = [];
    for i=1:length(pathwayID)
        regulatedPathwaysSorted = [regulatedPathwaysSorted, pathways(pathwayID(i)) ];
    end
    regulatedPathwaysSortedOverallTimepoints = regulatedPathwaysSorted';
    
    
    %%% reactions INCLUDING TRANSPORTERS
    % vector to save regulationvalue of each reaction 
    regulationReactions = ones(size(deviationMatrix,1),1)*NaN;
    % calculate overall regulation of reactions 
    for pathway = 1
        distance1 = abs(deviationMatrix(:,1) - deviationMatrix(:,2));
        distance2 = abs(deviationMatrix(:,2) - deviationMatrix(:,3));
        distance3 = abs(deviationMatrix(:,3) - deviationMatrix(:,4));
        overalldist = distance1 + distance2 + distance3;
        regulationReactions = overalldist;
        % rank them
        [sortvaluesReactions,rxnID] = sort(regulationReactions,'descend');
        regulatedReactionsSorted = [];
        longname = [];
        for i=1:length(rxnID)
            regulatedReactionsSorted = [regulatedReactionsSorted, rxnpathway(rxnID(i),2)];
            longname = [longname, rxnpathway(rxnID(i),3)];
        end
        regulatedReactionsSortedOveralltimepoints = [regulatedReactionsSorted',longname'];
    end
        
    %%% assign
    OV.overalltimepoints.regulatedPathwaysSortedOverallTimepoints = regulatedPathwaysSortedOverallTimepoints;
    OV.overalltimepoints.regulatedReactionsSortedOveralltimepoints = regulatedReactionsSortedOveralltimepoints;
    OV.overalltimepoints.sortvaluesPathway=sortvaluesPathway;
    OV.overalltimepoints.sortvaluesReactions=sortvaluesReactions;
    OV.overalltimepoints.pathwayID = pathwayID;
    OV.overalltimepoints.rxnID = rxnID;
    
    %%%%%%% timepoint 1to2
    
    %%% pathways
    % vector of overall regulation of each pathway in order to rank them
    regulationPathways = ones(size(pathways,1),1)*NaN;
    % calculate overall regulation of pathway     
    for pathway = 1:size(pathways,1)
        % assign pathwayspecific distancearray
        eval(strcat('x = pway.p',num2str(pathway),';'));
        % calculate distances for timepoints
        distance1 = sum(abs(x(:,1) - x(:,2)));
        regulationPathways(pathway) = distance1;
    end
    % rank them 
    [sortvaluesPathway,pathwayID] = sort(regulationPathways,'descend');
    regulatedPathwaysSorted = [];
    for i=1:length(pathwayID)
        regulatedPathwaysSorted = [regulatedPathwaysSorted, pathways(pathwayID(i)) ];
    end
    regulatedPathwaysSortedT12 = regulatedPathwaysSorted';
    % vector to save regulationvalue of each reaction 
    regulationReactions = ones(size(deviationMatrix,1),1)*NaN;
    % calculate overall regulation of reactions 
    for pathway = 1
        distance1 = abs(deviationMatrix(:,1) - deviationMatrix(:,2));
        overalldist = distance1;
        regulationReactions = overalldist;
        % rank them
        [sortvaluesReactions,rxnID] = sort(regulationReactions,'descend');
        regulatedReactionsSorted = [];
        longname = [];
        for i=1:length(rxnID)
            regulatedReactionsSorted = [regulatedReactionsSorted, rxnpathway(rxnID(i),2)];
            longname = [longname, rxnpathway(rxnID(i),3)];
        end
        regulatedReactionsSortedT12 = [regulatedReactionsSorted',longname'];
    end
    %%% assign 
    OV.t12.regulatedPathwaysSortedT12 = regulatedPathwaysSortedT12;
    OV.t12.regulatedReactionsSortedT12 = regulatedReactionsSortedT12;
    OV.t12.sortvaluesPathway=sortvaluesPathway;
    OV.t12.sortvaluesReactions=sortvaluesReactions; 
    OV.t12.pathwayID = pathwayID;
    OV.t12.rxnID = rxnID;
        
        
    
    %%%%%%% timepoint 2to3
    
    %%% pathways
    % vector of overall regulation of each pathway in order to rank them
    regulationPathways = ones(size(pathways,1),1)*NaN;
    % calculate overall regulation of pathway     
    for pathway = 1:size(pathways,1)
        % assign pathwayspecific distancearray
        eval(strcat('x = pway.p',num2str(pathway),';'));
        % calculate distances for timepoints
        distance2 = sum(abs(x(:,2) - x(:,3)));
        regulationPathways(pathway) = distance2;
    end
    % rank them 
    [sortvaluesPathway,pathwayID] = sort(regulationPathways,'descend');
    regulatedPathwaysSorted = [];
    for i=1:length(pathwayID)
        regulatedPathwaysSorted = [regulatedPathwaysSorted, pathways(pathwayID(i)) ];
    end
    regulatedPathwaysSortedT23 = regulatedPathwaysSorted';
    % vector to save regulationvalue of each reaction 
    regulationReactions = ones(size(deviationMatrix,1),1)*NaN;
    % calculate overall regulation of reactions 
    for pathway = 1
        distance2 = abs(deviationMatrix(:,2) - deviationMatrix(:,3));
        overalldist = distance2;
        regulationReactions = overalldist;
        % rank them
        [sortvaluesReactions,rxnID] = sort(regulationReactions,'descend');
        regulatedReactionsSorted = [];
        longname = [];
        for i=1:length(rxnID)
            regulatedReactionsSorted = [regulatedReactionsSorted, rxnpathway(rxnID(i),2)];
            longname = [longname, rxnpathway(rxnID(i),3)];
        end
        regulatedReactionsSortedT23 = [regulatedReactionsSorted',longname'];
    end
    %%% assign 
    OV.t23.regulatedPathwaysSortedT23 = regulatedPathwaysSortedT23;
    OV.t23.regulatedReactionsSortedT23 = regulatedReactionsSortedT23;    
    OV.t23.sortvaluesPathway=sortvaluesPathway;
    OV.t23.sortvaluesReactions=sortvaluesReactions;  
    OV.t23.pathwayID = pathwayID;
    OV.t23.rxnID = rxnID;
    
    
    %%%%%%% timepoint 3to4
    
    %%% pathways
    % vector of overall regulation of each pathway in order to rank them
    regulationPathways = ones(size(pathways,1),1)*NaN;
    % calculate overall regulation of pathway     
    for pathway = 1:size(pathways,1)
        % assign pathwayspecific distancearray
        eval(strcat('x = pway.p',num2str(pathway),';'));
        % calculate distances for timepoints
        distance3 = sum(abs(x(:,3) - x(:,4)));
        regulationPathways(pathway) = distance3;
    end
    % rank them 
    [sortvaluesPathway,pathwayID] = sort(regulationPathways,'descend');
    regulatedPathwaysSorted = [];
    for i=1:length(pathwayID)
        regulatedPathwaysSorted = [regulatedPathwaysSorted, pathways(pathwayID(i)) ];
    end
    regulatedPathwaysSortedT34 = regulatedPathwaysSorted';
    % vector to save regulationvalue of each reaction 
    regulationReactions = ones(size(deviationMatrix,1),1)*NaN;
    % calculate overall regulation of reactions 
    for pathway = 1
        distance3 = abs(deviationMatrix(:,3) - deviationMatrix(:,4));
        overalldist = distance3;
        regulationReactions = overalldist;
        % rank them
        [sortvaluesReactions,rxnID] = sort(regulationReactions,'descend');
        regulatedReactionsSorted = [];
        longname = [];
        for i=1:length(rxnID)
            regulatedReactionsSorted = [regulatedReactionsSorted, rxnpathway(rxnID(i),2)];
            longname = [longname, rxnpathway(rxnID(i),3)];
        end
        regulatedReactionsSortedT34 = [regulatedReactionsSorted',longname'];
    end
    %%% assign 
    OV.t34.regulatedPathwaysSortedT34 = regulatedPathwaysSortedT34;
    OV.t34.regulatedReactionsSortedT34 = regulatedReactionsSortedT34;     
    OV.t34.sortvaluesPathway=sortvaluesPathway;
    OV.t34.sortvaluesReactions=sortvaluesReactions;  
    OV.t34.pathwayID = pathwayID;
    OV.t34.rxnID = rxnID;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%% find overall top high and low regulated reactions per timepoint 
    %%% so basically just read in the already finished table 
    toprxnID = table2array(readtable('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/toprxnIDs.csv'));
    t0WTbiggerMTv = toprxnID(:,1);
    t0MTbiggerWTv = toprxnID(:,2);
    t1WTbiggerMTv = toprxnID(:,4);
    t1MTbiggerWTv = toprxnID(:,5);
    t2WTbiggerMTv = toprxnID(:,7);
    t2MTbiggerWTv = toprxnID(:,8);
    t3WTbiggerMTv = toprxnID(:,10);
    t3MTbiggerWTv = toprxnID(:,11);
    t0WTbiggerMT = [];
    t0MTbiggerWT = [];
    t1WTbiggerMT = [];
    t1MTbiggerWT = [];
    t2WTbiggerMT = [];
    t2MTbiggerWT = [];
    t3WTbiggerMT = [];
    t3MTbiggerWT = [];
    %%%%
    for i=1:size(toprxnID,1)
        %test = [test;rxnpathway(t0WTbiggerMT(i),1:4)];
        t0WTbiggerMT = [t0WTbiggerMT;rxnpathway(t0WTbiggerMTv(i),1:4)];
        t0MTbiggerWT = [t0MTbiggerWT;rxnpathway(t0MTbiggerWTv(i),1:4)];
        t1WTbiggerMT = [t1WTbiggerMT;rxnpathway(t1WTbiggerMTv(i),1:4)];
        t1MTbiggerWT = [t1MTbiggerWT;rxnpathway(t1MTbiggerWTv(i),1:4)];
        t2WTbiggerMT = [t2WTbiggerMT;rxnpathway(t2WTbiggerMTv(i),1:4)];
        t2MTbiggerWT = [t2MTbiggerWT;rxnpathway(t2MTbiggerWTv(i),1:4)];
        t3WTbiggerMT = [t3WTbiggerMT;rxnpathway(t3WTbiggerMTv(i),1:4)];
        t3MTbiggerWT = [t3MTbiggerWT;rxnpathway(t3MTbiggerWTv(i),1:4)];
    end
    OV.t1.t0WTbiggerMT = t0WTbiggerMT;
    OV.t1.t0MTbiggerWT = t0MTbiggerWT;
    OV.t2.t1WTbiggerMT = t1WTbiggerMT;
    OV.t2.t1MTbiggerWT = t1MTbiggerWT;
    OV.t3.t2WTbiggerMT = t2WTbiggerMT;
    OV.t3.t2MTbiggerWT = t2MTbiggerWT;
    OV.t4.t3WTbiggerMT = t3WTbiggerMT;
    OV.t4.t3MTbiggerWT = t3MTbiggerWT;
    
    
    
    
    %%% find rxn pathway with highest overall distances per timepoint 
    rankert0 = [];
    rankert1 = [];
    rankert2 = [];
    rankert3 = [];
    % assign sum of distances per timepoint per pathway 
    for pathway = 1:size(foundpathway,2)
        eval(strcat('x = pway.p',num2str(pathway),';'));
        rankert0(pathway) = sum(abs(x(:,1)));
        rankert1(pathway) = sum(abs(x(:,2)));
        rankert2(pathway) = sum(abs(x(:,3)));
        rankert3(pathway) = sum(abs(x(:,4)));
    end        
    % rank the pathwayregulation value per timepoint 
    [rankert0sorted,rankert0ID] = sort(rankert0,'descend');
    [rankert1sorted,rankert1ID] = sort(rankert1,'descend');
    [rankert2sorted,rankert2ID] = sort(rankert2,'descend');
    [rankert3sorted,rankert3ID] = sort(rankert3,'descend');
    % 
    % t0
    regulatedPathwaysT0 = [];
    for i=1:length(rankert0ID)
        regulatedPathwaysT0 = [regulatedPathwaysT0, pathways(rankert0ID(i)) ];
    end
    regulatedPathwaysSortedT0 = regulatedPathwaysT0';
    % t1
    regulatedPathwaysT1 = [];
    for i=1:length(rankert1ID)
        regulatedPathwaysT1 = [regulatedPathwaysT1, pathways(rankert1ID(i)) ];
    end
    regulatedPathwaysSortedT1 = regulatedPathwaysT1';
    % t2
    regulatedPathwaysT2 = [];
    for i=1:length(rankert2ID)
        regulatedPathwaysT2 = [regulatedPathwaysT2, pathways(rankert2ID(i)) ];
    end
    regulatedPathwaysSortedT2 = regulatedPathwaysT2';
    % t3
    regulatedPathwaysT3 = [];
    for i=1:length(rankert3ID)
        regulatedPathwaysT3 = [regulatedPathwaysT3, pathways(rankert3ID(i)) ];
    end
    regulatedPathwaysSortedT3 = regulatedPathwaysT3';
    % assign
    OV.t1.regulatedPathwaysSortedT0 = regulatedPathwaysSortedT0;
    OV.t2.regulatedPathwaysSortedT1 = regulatedPathwaysSortedT1;
    OV.t3.regulatedPathwaysSortedT2 = regulatedPathwaysSortedT2;
    OV.t4.regulatedPathwaysSortedT3 = regulatedPathwaysSortedT3;
    OV.t1.rankert0sorted = rankert0sorted;
    OV.t2.rankert1sorted = rankert1sorted;
    OV.t3.rankert2sorted = rankert2sorted;
    OV.t4.rankert3sorted = rankert3sorted;
    
    
    
    

end
   