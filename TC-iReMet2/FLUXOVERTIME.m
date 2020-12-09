%%% read in tables of fluxdeviations, reactions and their
%%% pathwayinformation, etc.
%%%%%%%%%%%%
%%% INIT %%%
%%%%%%%%%%%%

%%% read in tables
% read in rxnlist with corresponding pathways 
rxnpathway = readtable('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Data/RxnsList_ID_shortname_pathway.csv');
rxnpathway = table2cell(rxnpathway);
% read in deviationmatrix of fluxes between WTMT
deviationMatrix = dlmread('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/day/without_constCofactorRatios/with_slacks/with_slack_minimization/1e8/rep1/deviationMatrix_-_norm_to_Col0.dat','\t');
% norm deviationMatrix to dayspecific max in order to scale from -1 to 1
deviationMatrixNormed = ones(size(deviationMatrix))*NaN;
for col=1:size(deviationMatrix,2)
    deviationMatrixNormed(:,col) = deviationMatrix(:,col)/max(abs(deviationMatrix(:,col)));
end
% extract pathways and exclude empty cells 
pathways = rxnpathway(:,5);
pathways = pathways(~cellfun('isempty',pathways));
pathways = unique(pathways);
% create array where 1 corresponds to present reaction in pathway(column) 
foundpathway = ones(size(rxnpathway,1),size(pathways,1))*NaN;
for pathway = 1:size(pathways,1)
    foundpathway(:,pathway)=contains(rxnpathway(:,4),pathways{pathway});
end
% create table of fluxdeviation and reactionIDs for plotting 
for pathway = 1:size(pathways,1)
    eval(strcat ('pway.p' ,num2str(pathway), ' = [deviationMatrix(find(foundpathway(:,' ,num2str(pathway), ')),:), find(foundpathway(:,' ,num2str(pathway), '))];'));
end
% create table of fluxdeviations and reactionids of normed 
for pathway = 1:size(pathways,1)
    eval(strcat ('pwayNormed.p' ,num2str(pathway), ' = [deviationMatrixNormed(find(foundpathway(:,' ,num2str(pathway), ')),:), find(foundpathway(:,' ,num2str(pathway), '))];'));
end
% save on AI (allinfo)
AI.rxnpathway = rxnpathway;
AI.deviationMatrix = deviationMatrix;
AI.deviationMatrixNormed = deviationMatrixNormed;
AI.pathways = pathways;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATION OVER ALL TIMEPOINTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% pathways
% vector of overall regulation of each pathway in order to rank them
regulationValuePathways = ones(size(pathways,1),1)*NaN;
% calculate overall regulation of each pathway with manhatten distance     
for pathway = 1:size(pathways,1)
    % assign pathwayspecific distancearray
    eval(strcat('x = pway.p',num2str(pathway),';'));
    % calculate distances for timepoints
    distance1 = sum(abs(x(:,1) - x(:,2)));
    distance2 = sum(abs(x(:,2) - x(:,3)));
    distance3 = sum(abs(x(:,3) - x(:,4)));
    overalldist = distance1 + distance2 + distance3;
    % assign pathwayspecific overall regulation value
    regulationValuePathways(pathway) = overalldist/size(x,1); % needed to not skew it towards large amount of reactions(of course their net overall is bigger)
end
% sort them according to their overall regulation value 
[sortedValuesPathway,pathwayID] = sort(regulationValuePathways,'descend');
regulatedPathwaysSorted = [];
for id=1:length(pathwayID)
    regulatedPathwaysSorted = [regulatedPathwaysSorted, pathways(pathwayID(id)) ];
end
regulatedPathwaysSorted_OverAllTimepoints = regulatedPathwaysSorted';
% assign to AllInformation struc
AI.t_overall.pathways.sortedValuesPathway = sortedValuesPathway;
AI.t_overall.pathways.sortedPathwayID = pathwayID;
AI.t_overall.pathways.regulatedPathwaysSorted_OverAllTimepoints = regulatedPathwaysSorted_OverAllTimepoints;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% reactions
% calculate overall regulation overall timepoint for each reaction based on manhatten distance 
for pathway = 1
    distance1 = abs(deviationMatrix(:,1) - deviationMatrix(:,2));
    distance2 = abs(deviationMatrix(:,2) - deviationMatrix(:,3));
    distance3 = abs(deviationMatrix(:,3) - deviationMatrix(:,4));
    overalldist = distance1 + distance2 + distance3;
    regulationValueReactions = overalldist;
    % rank them
    [sortedvaluesReactions,rxnID] = sort(regulationValueReactions,'descend');
    regulatedReactionsSorted = [];
    longname = [];
    for i=1:length(rxnID)
        regulatedReactionsSorted = [regulatedReactionsSorted, rxnpathway(rxnID(i),2)];
        longname = [longname, rxnpathway(rxnID(i),3)];
    end
    regulatedReactionsSorted_OverAllTimepoints = [regulatedReactionsSorted',longname'];
end
% assign to AllInformation struc
AI.t_overall.reactions.sortedValuesReactions = sortedvaluesReactions;
AI.t_overall.reactions.sortedReactionID = rxnID;
AI.t_overall.reactions.regulatedReactionsSorted_OverAllTimepoints = regulatedReactionsSorted_OverAllTimepoints;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATION FOR EACH TIMEPOINTSTEP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% pathways %%%

%%% t0 to t1 (0d to 1d)
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
regulatedPathwaysSorted_T01 = regulatedPathwaysSorted';
% assign to AllInformation struc
AI.t_step_t01.pathways.sortedValuesPathway = sortvaluesPathway;
AI.t_step_t01.pathways.sortedPathwayID = pathwayID;
AI.t_step_t01.pathways.regulatedPathwaysSorted_t01 = regulatedPathwaysSorted_T01;

%%% t1 to t2 (1d to 7d)
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
regulatedPathwaysSorted_T12 = regulatedPathwaysSorted';
% assign to AllInformation struc
AI.t_step_t12.pathways.sortedValuesPathway = sortvaluesPathway;
AI.t_step_t12.pathways.sortedPathwayID = pathwayID;
AI.t_step_t12.pathways.regulatedPathwaysSorted_t12 = regulatedPathwaysSorted_T12;

%%% t2 to t3 (7d to 21d)
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
regulatedPathwaysSorted_T23 = regulatedPathwaysSorted';
% assign to AllInformation struc
AI.t_step_t23.pathways.sortedValuesPathway = sortvaluesPathway;
AI.t_step_t23.pathways.sortedPathwayID = pathwayID;
AI.t_step_t23.pathways.regulatedPathwaysSorted_t23 = regulatedPathwaysSorted_T23;


%%% reactions %%%

%%% t0 to t1 (0d to 1d)
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
    regulatedReactionsSorted_t01 = [regulatedReactionsSorted',longname'];
end
% assign to AllInformation struc
AI.t_step_t01.reactions.sortedValuesReactions = sortvaluesReactions;
AI.t_step_t01.reactions.sortedReactionsID = rxnID;
AI.t_step_t01.reactions.regulatedReactionsSorted_t01 = regulatedReactionsSorted_t01;

%%% t1 to t2 (1d to 7d)
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
    regulatedReactionsSorted_t12 = [regulatedReactionsSorted',longname'];
end
% assign to AllInformation struc
AI.t_step_t12.reactions.sortedValuesReactions = sortvaluesReactions;
AI.t_step_t12.reactions.sortedReactionsID = rxnID;
AI.t_step_t12.reactions.regulatedReactionsSorted_t12 = regulatedReactionsSorted_t12;

%%% t2 to t3 (7d to 21d)
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
    regulatedReactionsSorted_t23 = [regulatedReactionsSorted',longname'];
end
% assign to AllInformation struc
AI.t_step_t23.reactions.sortedValuesReactions = sortvaluesReactions;
AI.t_step_t23.reactions.sortedReactionsID = rxnID;
AI.t_step_t23.reactions.regulatedReactionsSorted_t23 = regulatedReactionsSorted_t23;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATION FOR EACH TIMEPOINT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% pathways %%%


%%% find rxn pathway with highest overall distances per timepoint 
rankert0 = [];
rankert1 = [];
rankert2 = [];
rankert3 = [];
% assign regulationsvalue per pathway 
for pathway = 1:size(foundpathway,2)
    eval(strcat('x = pway.p',num2str(pathway),';'));
    rankert0(pathway) = sum(abs(x(:,1)))/size(x,1);
    rankert1(pathway) = sum(abs(x(:,2)))/size(x,1);
    rankert2(pathway) = sum(abs(x(:,3)))/size(x,1);
    rankert3(pathway) = sum(abs(x(:,4)))/size(x,1);
end        
% rank the pathwayregulation value per timepoint 
[rankert0sorted,rankert0ID] = sort(rankert0,'descend');
[rankert1sorted,rankert1ID] = sort(rankert1,'descend');
[rankert2sorted,rankert2ID] = sort(rankert2,'descend');
[rankert3sorted,rankert3ID] = sort(rankert3,'descend');
%%%
regulatedPathwaysT0 = [];
regulatedPathwaysT1 = [];
regulatedPathwaysT2 = [];
regulatedPathwaysT3 = [];
for i=1:size(foundpathway,2)
    regulatedPathwaysT0 = [regulatedPathwaysT0, pathways(rankert0ID(i)) ];
    regulatedPathwaysT1 = [regulatedPathwaysT1, pathways(rankert1ID(i)) ];
    regulatedPathwaysT2 = [regulatedPathwaysT2, pathways(rankert2ID(i)) ];
    regulatedPathwaysT3 = [regulatedPathwaysT3, pathways(rankert3ID(i)) ];
end
regulatedPathwaysSortedT0 = regulatedPathwaysT0';
regulatedPathwaysSortedT1 = regulatedPathwaysT1';
regulatedPathwaysSortedT2 = regulatedPathwaysT2';
regulatedPathwaysSortedT3 = regulatedPathwaysT3';
%%% assign to big AI struc
AI.t_t0.pathways.sortedValuesPathway = rankert0sorted;
AI.t_t0.pathways.sortedPathwayID = rankert0ID';
AI.t_t0.pathways.regulatedPathwaysSorted_t0 = regulatedPathwaysSortedT0;
%
AI.t_t1.pathways.sortedValuesPathway = rankert1sorted;
AI.t_t1.pathways.sortedPathwayID = rankert1ID';
AI.t_t1.pathways.regulatedPathwaysSorted_t1 = regulatedPathwaysSortedT1;
%
AI.t_t2.pathways.sortedValuesPathway = rankert2sorted;
AI.t_t2.pathways.sortedPathwayID = rankert2ID';
AI.t_t2.pathways.regulatedPathwaysSorted_t2 = regulatedPathwaysSortedT2;
%
AI.t_t3.pathways.sortedValuesPathway = rankert3sorted;
AI.t_t3.pathways.sortedPathwayID = rankert3ID';
AI.t_t3.pathways.regulatedPathwaysSorted_t3 = regulatedPathwaysSortedT3;




%%% reactions %%%

% sort deviations per timepoint and find idx of reactions most regulated
% t0 
[WT_bigger_MT_t0, WT_bigger_MT_t0_ID] = sort(deviationMatrix(:,1),'descend');
[MT_bigger_WT_t0, MT_bigger_WT_t0_ID] = sort(deviationMatrix(:,1));
% t1 
[WT_bigger_MT_t1, WT_bigger_MT_t1_ID] = sort(deviationMatrix(:,2),'descend');
[MT_bigger_WT_t1, MT_bigger_WT_t1_ID] = sort(deviationMatrix(:,2));
% t2
[WT_bigger_MT_t2, WT_bigger_MT_t2_ID] = sort(deviationMatrix(:,3),'descend');
[MT_bigger_WT_t2, MT_bigger_WT_t2_ID] = sort(deviationMatrix(:,3));
% t3 
[WT_bigger_MT_t3, WT_bigger_MT_t3_ID] = sort(deviationMatrix(:,4),'descend');
[MT_bigger_WT_t3, MT_bigger_WT_t3_ID] = sort(deviationMatrix(:,4));
% make tables of top regulated names
% t0 WTMT
regulatedReactionsSorted = [];
longname = [];
for i=1:length(WT_bigger_MT_t0_ID)
    regulatedReactionsSorted = [regulatedReactionsSorted, rxnpathway(WT_bigger_MT_t0_ID(i),2)];
    longname = [longname, rxnpathway(WT_bigger_MT_t0_ID(i),3)];
end
regulatedReactionsSorted_WTbiggerMT = [regulatedReactionsSorted',longname'];
% t0 MTWT
regulatedReactionsSorted = [];
longname = [];
for i=1:length(WT_bigger_MT_t0_ID)
    regulatedReactionsSorted = [regulatedReactionsSorted, rxnpathway(MT_bigger_WT_t0_ID(i),2)];
    longname = [longname, rxnpathway(MT_bigger_WT_t0_ID(i),3)];
end
regulatedReactionsSorted_MTbiggerWT = [regulatedReactionsSorted',longname'];
%%% assign
AI.t_t0.reactions.sortedValuesReaction_WTbiggerMT = WT_bigger_MT_t0;
AI.t_t0.reactions.sortedValuesReaction_MTbiggerWT = MT_bigger_WT_t0;
AI.t_t0.reactions.sortedReactionID_WTbiggerMT = WT_bigger_MT_t0_ID;
AI.t_t0.reactions.sortedReactionID_MTbiggerWT = MT_bigger_WT_t0_ID;
AI.t_t0.reactions.regulatedReactionsSorted_WTbiggerMT = regulatedReactionsSorted_WTbiggerMT';
AI.t_t0.reactions.regulatedReactionsSorted_MTbiggerWT = regulatedReactionsSorted_MTbiggerWT';

% t1 WTMT
regulatedReactionsSorted = [];
longname = [];
for i=1:length(WT_bigger_MT_t0_ID)
    regulatedReactionsSorted = [regulatedReactionsSorted, rxnpathway(WT_bigger_MT_t1_ID(i),2)];
    longname = [longname, rxnpathway(WT_bigger_MT_t1_ID(i),3)];
end
regulatedReactionsSorted_WTbiggerMT = [regulatedReactionsSorted',longname'];
% t1 MTWT
regulatedReactionsSorted = [];
longname = [];
for i=1:length(WT_bigger_MT_t0_ID)
    regulatedReactionsSorted = [regulatedReactionsSorted, rxnpathway(MT_bigger_WT_t1_ID(i),2)];
    longname = [longname, rxnpathway(MT_bigger_WT_t1_ID(i),3)];
end
regulatedReactionsSorted_MTbiggerWT = [regulatedReactionsSorted',longname'];
%%% assign
AI.t_t1.reactions.sortedValuesReaction_WTbiggerMT = WT_bigger_MT_t1;
AI.t_t1.reactions.sortedValuesReaction_MTbiggerWT = MT_bigger_WT_t1;
AI.t_t1.reactions.sortedReactionID_WTbiggerMT = WT_bigger_MT_t1_ID;
AI.t_t1.reactions.sortedReactionID_MTbiggerWT = MT_bigger_WT_t1_ID;
AI.t_t1.reactions.regulatedReactionsSorted_WTbiggerMT = regulatedReactionsSorted_WTbiggerMT';
AI.t_t1.reactions.regulatedReactionsSorted_MTbiggerWT = regulatedReactionsSorted_MTbiggerWT';

% t2 WTMT
regulatedReactionsSorted = [];
longname = [];
for i=1:length(WT_bigger_MT_t0_ID)
    regulatedReactionsSorted = [regulatedReactionsSorted, rxnpathway(WT_bigger_MT_t2_ID(i),2)];
    longname = [longname, rxnpathway(WT_bigger_MT_t2_ID(i),3)];
end
regulatedReactionsSorted_WTbiggerMT = [regulatedReactionsSorted',longname'];
% t0 MTWT
regulatedReactionsSorted = [];
longname = [];
for i=1:length(WT_bigger_MT_t0_ID)
    regulatedReactionsSorted = [regulatedReactionsSorted, rxnpathway(MT_bigger_WT_t2_ID(i),2)];
    longname = [longname, rxnpathway(MT_bigger_WT_t2_ID(i),3)];
end
regulatedReactionsSorted_MTbiggerWT = [regulatedReactionsSorted',longname'];
%%% assign
AI.t_t2.reactions.sortedValuesReaction_WTbiggerMT = WT_bigger_MT_t2;
AI.t_t2.reactions.sortedValuesReaction_MTbiggerWT = MT_bigger_WT_t2;
AI.t_t2.reactions.sortedReactionID_WTbiggerMT = WT_bigger_MT_t2_ID;
AI.t_t2.reactions.sortedReactionID_MTbiggerWT = MT_bigger_WT_t2_ID;
AI.t_t2.reactions.regulatedReactionsSorted_WTbiggerMT = regulatedReactionsSorted_WTbiggerMT';
AI.t_t2.reactions.regulatedReactionsSorted_MTbiggerWT = regulatedReactionsSorted_MTbiggerWT';

% t3 WTMT
regulatedReactionsSorted = [];
longname = [];
for i=1:length(WT_bigger_MT_t0_ID)
    regulatedReactionsSorted = [regulatedReactionsSorted, rxnpathway(WT_bigger_MT_t3_ID(i),2)];
    longname = [longname, rxnpathway(WT_bigger_MT_t3_ID(i),3)];
end
regulatedReactionsSorted_WTbiggerMT = [regulatedReactionsSorted',longname'];
% t3 MTWT
regulatedReactionsSorted = [];
longname = [];
for i=1:length(WT_bigger_MT_t0_ID)
    regulatedReactionsSorted = [regulatedReactionsSorted, rxnpathway(MT_bigger_WT_t3_ID(i),2)];
    longname = [longname, rxnpathway(MT_bigger_WT_t3_ID(i),3)];
end
regulatedReactionsSorted_MTbiggerWT = [regulatedReactionsSorted',longname'];
%%% assign
AI.t_t3.reactions.sortedValuesReaction_WTbiggerMT = WT_bigger_MT_t3;
AI.t_t3.reactions.sortedValuesReaction_MTbiggerWT = MT_bigger_WT_t3;
AI.t_t3.reactions.sortedReactionID_WTbiggerMT = WT_bigger_MT_t3_ID;
AI.t_t3.reactions.sortedReactionID_MTbiggerWT = MT_bigger_WT_t3_ID;
AI.t_t3.reactions.regulatedReactionsSorted_WTbiggerMT = regulatedReactionsSorted_WTbiggerMT';
AI.t_t3.reactions.regulatedReactionsSorted_WTbiggerMT = regulatedReactionsSorted_MTbiggerWT';


%%% without transporter
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
AI.t_t0.reactions.top15_noTransporter_WTbiggerMT = t0WTbiggerMT;
AI.t_t0.reactions.top15_noTransporter_MTbiggerWT = t0MTbiggerWT;
AI.t_t1.reactions.top15_noTransporter_WTbiggerMT = t1WTbiggerMT;
AI.t_t1.reactions.top15_noTransporter_MTbiggerWT = t1MTbiggerWT;
AI.t_t2.reactions.top15_noTransporter_WTbiggerMT = t2WTbiggerMT;
AI.t_t2.reactions.top15_noTransporter_MTbiggerWT = t2MTbiggerWT;
AI.t_t3.reactions.top15_noTransporter_WTbiggerMT = t3WTbiggerMT;
AI.t_t3.reactions.top15_noTransporter_MTbiggerWT = t3MTbiggerWT;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%    PLOTTING    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% for plotting normalized distances are required for this use the
%%% normedDevmatrix to plot graphs for them to be correct and compareable

%%% create vertical barplots of flux deviations
for pathway = 1:size(foundpathway,2)
    eval(strcat('x = pwayNormed.p',num2str(pathway),';'));
    bar(x(:,1:4));
    title(strcat('Reaction pathway: ',' ',pathways{pathway}));
    set(gca,'XTickLabel',(rxnpathway(x(:,5),3)'));
    xtickangle(45);
    ylim([-1 1]);
    %set(gca,'YScale','log');
    ylabel('Distance in [ mol^2 / (gDW d)^2]');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/pathway_normedToTimepointSpecificAbsMaxValue/pathway',num2str(pathway),'_',regexprep(pathways{pathway},'/',''),'.png'));
end

for pathway = 1:size(foundpathway,2)
    eval(strcat('x = pwayNormed.p',num2str(pathway),';'));
    bar(x(:,1:4));
    title(strcat('Reaction pathway: ',' ',pathways{pathway}));
    set(gca,'XTickLabel',(rxnpathway(x(:,5),2)'));
    xtickangle(45);
    ylim([-1 1]);
    %set(gca,'YScale','log');
    ylabel('Distance in [ mol^2 / (gDW d)^2]');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/pathway_normedToTimepointSpecificAbsMaxValue_shortNames/pathway',num2str(pathway),'_',regexprep(pathways{pathway},'/',''),'.png'));
end



%%%% try normal values n stuff
for pathway = 1:size(foundpathway,2)
    eval(strcat('x = pway.p',num2str(pathway),';'));
    bar(x(:,1:4));
    title(strcat('Reaction pathway: ',' ',pathways{pathway}));
    set(gca,'XTickLabel',(rxnpathway(x(:,5),3)'));
    xtickangle(45);
    ylim([-30 30]);
    %set(gca,'YScale','log');
    ylabel('Distance in [ mol^2 / (gDW d)^2]');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/pathway_scaleInValueRange30/pathway',num2str(pathway),'_',regexprep(pathways{pathway},'/',''),'.png'));
end

%%%% try log scale n stuff
for pathway = 1:size(foundpathway,2)
    eval(strcat('x = pway.p',num2str(pathway),';'));
    bar(abs(x(:,1:4)));
    b = bar(x(:,1:4),'FaceColor','flat');
    title(strcat('Reaction pathway: ',' ',pathways{pathway}));
    set(gca,'XTickLabel',(rxnpathway(x(:,5),3)'));
    xtickangle(45);
    %ylim([-30 30]);
    %set(gca,'YScale','log');
    ylabel('Distance in [ mol^2 / (gDW d)^2]');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/pathway_scaleInValueRange30/pathway',num2str(pathway),'_',regexprep(pathways{pathway},'/',''),'.png'));
end



pway.p1
    





figure
x1 = 1:size(x,1);
y1 = x(:,1);
bar(x1,y1)
ax1 = gca; % current axes
ax1.XColor = 'r';
ax1.YColor = 'r';

ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');

x2 = 1:size(x,1);
y2 = x(:,1);
bar(x2,y2,'Parent',ax2)




figure
l1=[1:10];
l2=[1:10];
a=[rand(10,1)*100 zeros(10,1)  ];
b=[zeros(10,1)    rand(10,1)*5 ];
[AX,H1,H2] = plotyy(l1,a, l2,b, 'bar', 'bar');
set(H1,'FaceColor','r') % a
set(H2,'FaceColor','b') % b
hold on
[AX,H1,H2] = plotyy([1:10],a, [1:10],b, 'bar', 'bar');
set(H1,'FaceColor','g') % a
set(H2,'FaceColor','y') % b


n = 1:10;
x = [1:10];
y = [1:10];
a = [1 NaN 3 4];
b = abs([-1 -2 -3 -4]);
c = [NaN 2 NaN NaN; NaN NaN NaN NaN];

 

yyaxis left
ylim([-10 10])

bar([a;b]);

ax1 = gca;
set(ax,'YScale','log')
ylim('log')

yyaxis right
ylim([-10 10])

set(gca,'YScale','log')
ax = gca;
ax.YDir = 'reverse'

bar(c)


n = 1:10;
x = [1:10];
y = [1:10];
a = [1 NaN 3 4];
b = abs([-1 -2 -3 -4]);
c = [NaN 2 NaN NaN; NaN NaN NaN NaN];




yyaxis left
bar([a;b]);
set(gca,'YScale','log')





yyaxis right
bar(c)
ax2 = gca;
set(ax2,'YScale','log')
gca.YDir ='reverse';





