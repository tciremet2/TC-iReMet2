function [] = Eval_plotTopRegulatedPathways

	clear;
	%%% read in data and set important variables
	kNumbRxns = 549;
	rxnpathway = readtable('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Data/RxnsList_ID_shortname_pathway.csv');
    % read in deviationmatrix of fluxes between WTMT
    deviationMatrix = dlmread('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/day/without_constCofactorRatios/with_slacks/with_slack_minimization/1e8/rep1/deviationMatrix_-_norm_to_Col0.dat','\t');
    deviationMatrix(abs(deviationMatrix)<1e-4) = 0;
	load artificalRxns.mat
	
    %%% 
    % convert to cell for easier handling 
    rxnpathway = table2cell(rxnpathway);
    % extract pathways and exclude empty cells 
    pathways = rxnpathway(:,5);
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
        eval(strcat ('pway.withTnA.p' ,num2str(pathway), '.M = [deviationMatrix(find(foundpathway(:,' ,num2str(pathway), ')),:), find(foundpathway(:,' ,num2str(pathway), '))];'));
        eval(strcat ('pway.withTnA.p' ,num2str(pathway), '.pwayname = pathways(', num2str(pathway), ');'));
	end
	
	%%%
    % create table of fluxdeviation and reactionIDs for plotting without t
	%deviationMatrix = deviationMatrix(setdiff(1:549,artificial_rxns_and_transporters)',:);
	foundpathway = ones(size(rxnpathway,1),size(pathways,1))*NaN;
	foundpathway = foundpathway(setdiff(1:549,artificial_rxns_and_transporters)',:);
	rxnpathway = rxnpathway(setdiff(1:549,artificial_rxns_and_transporters)',:);
    for pathway = 1:size(pathways,1)
        foundpathway(:,pathway)=contains(rxnpathway(:,4),pathways{pathway});
    end
    for pathway = 1:size(foundpathway,2)
        eval(strcat ('pway.withoutTnA.p' ,num2str(pathway), '.M = [deviationMatrix(find(foundpathway(:,' ,num2str(pathway), ')),:), find(foundpathway(:,' ,num2str(pathway), '))];'));
        eval(strcat ('pway.withoutTnA.p' ,num2str(pathway), '.pwayname = pathways(', num2str(pathway), ');'));
	end

	
	
	
	%%% revert to original 
	rxnpathway = readtable('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Data/RxnsList_ID_shortname_pathway.csv');
    % read in deviationmatrix of fluxes between WTMT
    deviationMatrix = dlmread('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/day/without_constCofactorRatios/with_slacks/with_slack_minimization/1e8/rep1/deviationMatrix_-_norm_to_Col0.dat','\t');
    deviationMatrix(abs(deviationMatrix)<1e-4) = 0;
	rxnpathway = table2cell(rxnpathway);
    pathways = rxnpathway(:,5);
    pathways = pathways(~cellfun('isempty',pathways));
    foundpathway = ones(size(rxnpathway,1),size(pathways,1))*NaN;
    for pathway = 1:size(pathways,1)
        foundpathway(:,pathway)=contains(rxnpathway(:,4),pathways{pathway});
    end
	
	
	%{
	for i=1:length(pathways)
		eval(strcat('disp(pway.p',num2str(i),'.pwayname)'))
	end
	%}

	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% FIND TOP REGULATED PATHWAYS ACCORDING TO MEASURE %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% with transportreactions %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%% MANHATTAN
	% create vector to assign ranking to 
	regulationPathways = ones(size(pathways,1),1)*NaN;
    % calculate overall regulation value of each pathway     
    for pathway = 1:size(pathways,1)
        % assign pathwayspecific distancearray
        eval(strcat('x = pway.withTnA.p',num2str(pathway),'.M;'));
		% normalize to row to access relative change
		x(:,1:4) = x(:,1:4)./max(abs(x(:,1:4)),[],2);
		x(isnan(x)) = 0;
        % calculate distances for timepoints
        distance1 = sum(abs(x(:,1) - x(:,2)));
        distance2 = sum(abs(x(:,2) - x(:,3)));
        distance3 = sum(abs(x(:,3) - x(:,4)));
        overalldist = distance1 + distance2 + distance3;
        % assign pathwayspecific overall regulation value
        regulationPathways(pathway) = overalldist/size(x,1);
    end
    % rank them 
    [sortvaluesPathway,pathwayID] = sort(regulationPathways,'descend');
    regulatedPathwaysSorted = [];
    for i=1:length(pathwayID)
        regulatedPathwaysSorted = [regulatedPathwaysSorted, pathways(pathwayID(i)) ];
    end
    regulatedPathwaysSortedOverallTimepoints = regulatedPathwaysSorted';	
	eval(strcat('pway.withTnA.Manhattan = regulatedPathwaysSortedOverallTimepoints;'));
	withTnA_man = [sortvaluesPathway, pathwayID];
	
	
	
	%%% SUM
	% create vector to assign ranking to 
	regulationPathways = ones(size(pathways,1),1)*NaN;
    % calculate overall regulation value of each pathway     
    for pathway = 1:size(pathways,1)
        % assign pathwayspecific distancearray
        eval(strcat('x = pway.withTnA.p',num2str(pathway),'.M;'));
		% normalize to row to access relative change
		x(:,1:4) = x(:,1:4)./max(abs(x(:,1:4)),[],2);
		x(isnan(x)) = 0;
        % calculate distances for timepoints
        distancevec = x(:,1:4);
		overalldist = sum(distancevec(:));        
        % assign pathwayspecific overall regulation value
        regulationPathways(pathway) = overalldist/size(x,1);
    end
    % rank them 
    [sortvaluesPathway,pathwayID] = sort(regulationPathways,'descend');
    regulatedPathwaysSorted = [];
    for i=1:length(pathwayID)
        regulatedPathwaysSorted = [regulatedPathwaysSorted, pathways(pathwayID(i)) ];
    end
    regulatedPathwaysSortedOverallTimepoints = regulatedPathwaysSorted';	
	eval(strcat('pway.withTnA.Sum = regulatedPathwaysSortedOverallTimepoints;'));
	withTnA_sum = [sortvaluesPathway, pathwayID];
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% without transportreactions %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%% MANHATTAN
	% create vector to assign ranking to 
	regulationPathways = ones(size(pathways,1),1)*NaN;
    % calculate overall regulation value of each pathway     
    for pathway = 1:size(pathways,1)
        % assign pathwayspecific distancearray
        eval(strcat('x = pway.withoutTnA.p',num2str(pathway),'.M;'));
		% normalize to row to access relative change
		x(:,1:4) = x(:,1:4)./max(abs(x(:,1:4)),[],2);
		x(isnan(x)) = 0;
        % calculate distances for timepoints
        distance1 = sum(abs(x(:,1) - x(:,2)));
        distance2 = sum(abs(x(:,2) - x(:,3)));
        distance3 = sum(abs(x(:,3) - x(:,4)));
        overalldist = distance1 + distance2 + distance3;
        % assign pathwayspecific overall regulation value
        regulationPathways(pathway) = overalldist/size(x,1);
    end
    % rank them 
    [sortvaluesPathway,pathwayID] = sort(regulationPathways,'descend');
    regulatedPathwaysSorted = [];
    for i=1:length(pathwayID)
        regulatedPathwaysSorted = [regulatedPathwaysSorted, pathways(pathwayID(i)) ];
    end
    regulatedPathwaysSortedOverallTimepoints = regulatedPathwaysSorted';	
	eval(strcat('pway.withoutTnA.Manhattan = regulatedPathwaysSortedOverallTimepoints;'));
	withoutTnA_man = [sortvaluesPathway, pathwayID];
	
	%%% SUM
	% create vector to assign ranking to 
	regulationPathways = ones(size(pathways,1),1)*NaN;
    % calculate overall regulation value of each pathway     
    for pathway = 1:size(pathways,1)
        % assign pathwayspecific distancearray
        eval(strcat('x = pway.withoutTnA.p',num2str(pathway),'.M;'));
		% normalize to row to access relative change
		x(:,1:4) = x(:,1:4)./max(abs(x(:,1:4)),[],2);
		x(isnan(x)) = 0;
        % calculate distances for timepoints
        distancevec = x(:,1:4);
		overalldist = sum(distancevec(:));        
        % assign pathwayspecific overall regulation value
        regulationPathways(pathway) = overalldist/size(x,1);
    end
    % rank them 
    [sortvaluesPathway,pathwayID] = sort(regulationPathways,'descend');
    regulatedPathwaysSorted = [];
    for i=1:length(pathwayID)
        regulatedPathwaysSorted = [regulatedPathwaysSorted, pathways(pathwayID(i)) ];
    end
    regulatedPathwaysSortedOverallTimepoints = regulatedPathwaysSorted';	
	eval(strcat('pway.withoutTnA.Sum = regulatedPathwaysSortedOverallTimepoints;'));
	withoutTnA_sum = [sortvaluesPathway, pathwayID];
	% create new naming 
	
	
	
	
	%%%%%%%%%%%%%%%%%%%%%
	%%% PLOT PATHWAYS %%%
	%%%%%%%%%%%%%%%%%%%%%
	
	
	%%%%%%%%%%%%%%%%%%%%%
	%%% plot with TnA %%%
	%%%%%%%%%%%%%%%%%%%%%
	
	%%% with TnA
	for pathway = 1:size(pathways,1)
		
		% extract pathway matrix 
		eval(strcat('x = pway.withTnA.p', num2str(pathway), '.M(:,1:4);'));
		eval(strcat('xb = pway.withTnA.p', num2str(pathway), '.M;'));
		x = sortrows(x,4,'descend');
		xb = sortrows(xb,4,'descend');
		Ndecimals = 2; 
		f = 10.^Ndecimals ;
		for i = 1:size(x,1)
			for j = 1:size(x,2)
				x(i,j) = round(f*x(i,j))/f;
			end
		end
		% cell trans
		n = num2cell(xb);
		% exchange ids with actual reaction names 
		for i = 1:size(n,1)
			id = n{i,5};
			n{i,5} = rxnpathway{id,2};
		end
		% convert to table for heatmap 
		n = cell2table(n,'VariableNames',{'day_0','day_1','day_7','day_21','nmbr'});
		% reorganize the table 
		c =1; % counter
		crxnnames = strings(1,length(x(:)));
		ctimepoint = strings(1,length(x(:)));
		cvalue = ones(1,length(x(:)))*NaN;
		for row = 1:size(x,1)
			for col = 1:size(x,2)			
				crxnnames(c) = n{row,5};
				ctimepoint(c) = string(n.Properties.VariableNames{col});
				cvalue(c) = x(row,col);
				% update c			
				c = c + 1;		
			end
		end
		T = table(crxnnames(:),ctimepoint(:),cvalue(:),'VariableNames',{'reaction','timepoint','fluxdiff'});
		h = heatmap(T,'timepoint','reaction','ColorVariable','fluxdiff','FontSize',13);
		oi = strcat('h.Title = ''pathway:',{' '} ,pathways(pathway) ,''';');
		eval(oi{1})
		%h.XLabel = 'timepoints';
		h.XLabel = '';
		h.YLabel = '';
		%h.FontSize = 12;
		h.ColorLimits = [-5 5];
		colormap(redblue())
		h.SourceTable.timepoint = categorical(h.SourceTable.timepoint);
		neworder = {'day_0','day_1','day_7','day_21'};
		h.SourceTable.timepoint = reordercats(h.SourceTable.timepoint,neworder);
		h.SourceTable.reaction = categorical(h.SourceTable.reaction);
		neworder2 = table2cell(n(:,end));
		h.SourceTable.reaction = reordercats(h.SourceTable.reaction,neworder2);
		%%% plot
		saveas(h,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/pathways/withTnA/pathway',num2str(pathway),'_',regexprep(pathways{pathway},'/',''),'.png'));

	end
			
		
	%%% without TnA
	for pathway = 1:size(pathways,1)
		
		% extract pathway matrix 
		eval(strcat('x = pway.withoutTnA.p', num2str(pathway), '.M(:,1:4);'));
		eval(strcat('xb = pway.withoutTnA.p', num2str(pathway), '.M;'));
		x = sortrows(x,4,'descend');
		xb = sortrows(xb,4,'descend');
		Ndecimals = 2; 
		f = 10.^Ndecimals ;
		for i = 1:size(x,1)
			for j = 1:size(x,2)
				x(i,j) = round(f*x(i,j))/f;
			end
		end
		if isempty(x)
			continue
		end		
		% cell trans
		n = num2cell(xb);
		% exchange ids with actual reaction names 
		for i = 1:size(n,1)
			id = n{i,5};
			n{i,5} = rxnpathway{id,2};
		end
		% convert to table for heatmap 
		n = cell2table(n,'VariableNames',{'day_0','day_1','day_7','day_21','nmbr'});
		% reorganize the table 
		c =1; % counter
		crxnnames = strings(1,length(x(:)));
		ctimepoint = strings(1,length(x(:)));
		cvalue = ones(1,length(x(:)))*NaN;
		for row = 1:size(x,1)
			for col = 1:size(x,2)			
				crxnnames(c) = n{row,5};
				ctimepoint(c) = string(n.Properties.VariableNames{col});
				cvalue(c) = x(row,col);
				% update c			
				c = c + 1;		
			end
		end
		T = table(crxnnames(:),ctimepoint(:),cvalue(:),'VariableNames',{'reaction','timepoint','fluxdiff'});
		h = heatmap(T,'timepoint','reaction','ColorVariable','fluxdiff','FontSize',13);
		oi = strcat('h.Title = ''pathway:',{' '} ,pathways(pathway) ,''';');
		eval(oi{1})
		%h.XLabel = 'timepoints';
		h.XLabel = '';
		h.YLabel = '';
		%h.FontSize = 12;
		h.ColorLimits = [-5 5];
		colormap(redblue())
		h.SourceTable.timepoint = categorical(h.SourceTable.timepoint);
		neworder = {'day_0','day_1','day_7','day_21'};
		h.SourceTable.timepoint = reordercats(h.SourceTable.timepoint,neworder);
		h.SourceTable.reaction = categorical(h.SourceTable.reaction);
		neworder2 = table2cell(n(:,end));
		h.SourceTable.reaction = reordercats(h.SourceTable.reaction,neworder2);
		%%% plot
		saveas(h,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/pathways/withoutTnA/pathway',num2str(pathway),'_',regexprep(pathways{pathway},'/',''),'.png'));

	end
	
	%{
	%%% plot overall deviation
    overalldev = [5e-02, 1.081e+01, 4.586e+01, 9.086e+01];
    bar((overalldev));
    %title('Distance over all timepoints: 0d-1d-7d-21d');
    set(gca,'XTickLabel',{'day 0','day 1','day 7','day 21'},'FontSize',27,'FontWeight','bold');
    set(gca,'YScale','log')
	ylim([1e-2 1e3])
    ylabel('Distance [ mol^2 / (gDW d)^2]','FontSize',30,'FontWeight','bold');
	xlabel('Timepoints','FontSize',30,'FontWeight','bold')
	%set(gca, 'XGrid','off');
	%set(gca, 'YGrid','off');
	xh = get(gca,'xlabel') % handle to the label object
	p = get(xh,'position') % get the current position property
	p(2) = 0.9*p(2) ;        % double the distance, 
                       % negative values put the label below the axis
	set(xh,'position',p)   % set the new position
	yh = get(gca,'ylabel') % handle to the label object
	p = get(yh,'position') % get the current position property
	p(1) = 1.2*p(1) ;        % double the distance, 
                       % negative values put the label below the axis
	set(yh,'position',p)   % set the new position
    set(gca, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
	text(1:length(overalldev),overalldev,num2str(overalldev'),'vert','bottom','horiz','center','Fontsize',40); 
    saveas(gca,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/Overall_Distance.png'));
	%}
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% CREATE RANKING INFORMATION TABLE AND SAVE IT %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%% init object
	ranking = [];
	ranking.withTnA_man = withTnA_man;
	ranking.withTnA_sum = withTnA_sum;
	ranking.withoutTnA_man = withoutTnA_man;
	ranking.withoutTnA_sum = withoutTnA_sum;
	ranking.withTnA_man = [num2cell(ranking.withTnA_man), pway.withTnA.Manhattan];
	ranking.withTnA_sum = [num2cell(ranking.withTnA_sum), pway.withTnA.Sum];
	ranking.withoutTnA_man = [num2cell(ranking.withoutTnA_man), pway.withoutTnA.Manhattan];
	ranking.withoutTnA_sum = [num2cell(ranking.withoutTnA_sum), pway.withoutTnA.Sum];
	% write to table 
	writetable(cell2table(ranking.withTnA_man),'/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/ranking_withTnA_man.csv')
	writetable(cell2table(ranking.withTnA_sum),'/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/ranking_withTnA_sum.csv')
	writetable(cell2table(ranking.withoutTnA_man),'/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/ranking_withoutTnA_man.csv')
	writetable(cell2table(ranking.withoutTnA_sum),'/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/ranking_withoutTnA_sum.csv')
	% to make new names according to rxn[compartment]
	% use Set_nameChanger 
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% EVAL TOP RXNS PER CLUSTER %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% read in cluster data
	clust = [];
	for i = 1:7	
		eval(strcat('clust.cluster',num2str(i),' = readtable(''/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/R_calcs/cluster', num2str(i), '.txt'');'));
		eval(strcat('clust.cluster',num2str(i),'=table2cell(clust.cluster',num2str(i),');'));
	end
	
	% reshape matrices to be heatmapable 		
	for cl = 1:7	
		% extract matrices\information
		eval(strcat('x = clust.cluster', num2str(cl), '(1:10,2:5);'));
		eval(strcat('xb = clust.cluster', num2str(cl), '(1:10,[2 3 4 5 7]);'));
		Ndecimals = 2; 
		f = 10.^Ndecimals ;
		for i = 1:size(x,1)
			for j = 1:size(x,2)
				x{i,j} = round(f*x{i,j})/f;
			end
		end
		if isempty(x)
			continue
		end		
		n = xb;
		for i = 1:size(n,1)
			id = n{i,5};
			n{i,5} = rxnpathway{id,2};
		end
		n = cell2table(n,'VariableNames',{'day_0','day_1','day_7','day_21','nmbr'});
		% reorganize the table 
		c =1; % counter
		crxnnames = strings(1,length(x(:)));
		ctimepoint = strings(1,length(x(:)));
		cvalue = ones(1,length(x(:)))*NaN;
		for row = 1:size(x,1)
			for col = 1:size(x,2)			
				crxnnames(c) = n{row,5};
				ctimepoint(c) = string(n.Properties.VariableNames{col});
				cvalue(c) = x{row,col};
				% update c			
				c = c + 1;		
			end
		end
		T = table(crxnnames(:),ctimepoint(:),cvalue(:),'VariableNames',{'reaction','timepoint','fluxdiff'});
		h = heatmap(T,'timepoint','reaction','ColorVariable','fluxdiff','FontSize',13);
		oi = strcat('h.Title = ''Cluster:',{' '} ,num2str(cl) ,''';');
		eval(oi{1})
		%h.XLabel = 'timepoints';
		h.XLabel = '';
		h.YLabel = '';
		%h.FontSize = 12;
		h.ColorLimits = [-5 5];
		colormap(redblue())
		h.SourceTable.timepoint = categorical(h.SourceTable.timepoint);
		neworder = {'day_0','day_1','day_7','day_21'};
		h.SourceTable.reaction = categorical(h.SourceTable.reaction);
		neworder2 = table2cell(n(:,end));
		h.SourceTable.timepoint = reordercats(h.SourceTable.timepoint,neworder);
		h.SourceTable.reaction = reordercats(h.SourceTable.reaction,neworder2);
		% save
		saveas(h,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/HEATMAPSOFKMEANSCLUSTERING/hmapCluster',num2str(cl),'.png'));
	end
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%	FIND EARLY TIMEPOINTS TOP REACTIONS %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%% top diff rxn 0d to 1d
	x = deviationMatrix;
	x = [x,[1:549]'];
	x = x(setdiff(1:549,artificial_rxns_and_transporters)',:);
	%x(:,1:4) = x(1:end,1:4)./max(abs(x(1:end,1:4)),[],2);
	%x(isnan(x))=0;
	x = [abs(x(:,1)-x(:,2)),x(:,5)];
	x = sortrows(x,'descend');
	x = [x,[1:320]'];
	n = num2cell(x);
	n = Set_nameChanger(n,rxnpathway);
	for i = 1:size(n,1)
			id = n{i,2};
			n{i,2} = rxnpathway{id,2};
			%n{i,3} = rxnpathway{id,3};
	end
	top0to1regRxns = n(1:15,:);	
	%%% find idx of top 10 for 0d and 1d
	x = deviationMatrix;
	x = [x,[1:549]'];
	x = x(:,[1 2 5]);
	x = x(setdiff(1:549,artificial_rxns_and_transporters)',:);
	[~,ID1]=sortrows(abs(x),1,'descend');
	[~,ID2]=sortrows(abs(x),2,'descend');
	x0d = x(ID1,:);
	top0d = num2cell(x0d(1:10,:));
	for i = 1:size(top0d,1)
			id = top0d{i,3};
			top0d{i,2} = rxnpathway{id,2};
			top0d{i,3} = rxnpathway{id,3};
	end
	x1d = x(ID2,:);
	x1d(:,[2 1 3]);
	top1d = num2cell(x1d(1:10,:));
	for i = 1:size(top1d,1)
			id = top1d{i,3};
			top1d{i,2} = rxnpathway{id,2};
			top1d{i,3} = rxnpathway{id,3};
	end
	% create cell
	early = [top0d(:,2),top1d(:,2),top0to1regRxns(1:10,3)];
	writetable(cell2table(early),'/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/earlytimepointstop.csv')

	
	
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% CREATE RANKING OF ALL TRANSPORT AND ARTIFICIAL REACTIONS %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%% rank artificial 
	x = deviationMatrix;
	x = [x,[1:549]'];
	%x = [x,[1:549]'];
	x = x(artificial_rxns_and_transporters,:);
	distance1 = (abs(x(:,1) - x(:,2)));
    distance2 = (abs(x(:,2) - x(:,3)));
    distance3 = (abs(x(:,3) - x(:,4)));
    overalldist = distance1 + distance2 + distance3;
	x = [x, overalldist];
	x = sortrows(x,6,'descend');
	x = [x,[1:size(x,1)]'];
	x = x(:,[1 2 3 4 5 7 6]);
	x = num2cell(x);
	%lad = Set_nameChanger(x(:,4:6),rxnpathway);
	for i = 1:size(x,1)
			id = x{i,5};
			x{i,5} = rxnpathway{id,2};
			x{i,6} = rxnpathway{id,3};
	end
	writetable(cell2table(x),'/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/topRxnTnA.csv')
	
	
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% U(p),D(own),N(o change) evaluation %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%% create udn cells %%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% udnclust: udncode, rxnid                 %%%
	%%% summary: clIDudncodes, udncode, amount   %%%
	%%% IdUdnClnmbr: RxnID, udncode, clIDudnCode %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	
	
	% read in cluster data
	clust = [];
	for i = 1:7	
		eval(strcat('clust.cluster',num2str(i),' = readtable(''/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/R_calcs/cluster', num2str(i), '.txt'');'));
		eval(strcat('clust.cluster',num2str(i),'=table2cell(clust.cluster',num2str(i),');'));
	end	

	% calculate and updowntable based on size of frame
	UDN = [];
	% 		
	for cl = 1:7	
		% extract matrices\information % 7rxnid, 2:5 values
		eval(strcat('x = clust.cluster', num2str(cl), '(:,[2 3 4 5 7]);'));		
		% calculate clusterspecific udn and save according to cluster		
		eval(strcat('udnclust',num2str(cl),'=cell(size(x,1),3);'));
		udnclust = cell(size(x,1),3);
		for row = 1:size(x,1)			
			udnclust{row,1} = x{row,2}-x{row,1};
			udnclust{row,2} = x{row,3}-x{row,2};
			udnclust{row,3} = x{row,4}-x{row,3};
		end
		% make changes that are basically 0 to 0
		udnclust = cell2mat(udnclust);
		udnclust(abs(udnclust)<1e-4) = 0;
		udnclust = num2cell(udnclust);
		% exchange with u d n for how the change behaves 
		for row = 1:size(udnclust,1)
			for col = 1:size(udnclust,2)				
				% check for posneg0value
				if udnclust{row,col}>0
					udnclust{row,col} = "U";
				elseif udnclust{row,col}<0
					udnclust{row,col} = "D";
				elseif isequal(udnclust{row,col},0)
					udnclust{row,col} = "N";
				end
				
			end
		end		
		% add id
		udnclust = [udnclust,x(:,end)];
		% save udn of that specific cluster 
		eval(strcat('UDN.udnCluster',num2str(cl),'.udnclust=udnclust;'));			
	
		%%% tighten frame to classes 
		udnclustNew = cell(size(udnclust,1),2);
		for row = 1:size(udnclust,1)
			udnclustNew{row,1} = strcat(udnclust{row,1},udnclust{row,2},udnclust{row,3});
		end
		udnclustNew(:,2)=udnclust(:,4);
		% save
		eval(strcat('UDN.udnCluster',num2str(cl),'.udnclust=udnclustNew;'));
		
			%%%
		%%% evaluate udn by 
		[reg,~,ic] = unique(table2array(cell2table(udnclustNew(:,1))));
		occurences = [reg,accumarray(ic,1)];
		gesamt = [	table2array(cell2table(udnclustNew(:,1))),ic];
		%gesamt = [ic,table2array(cell2table(udnclustNew))];
		key = unique(gesamt,'rows');
		% sort
		occurences = sortrows(occurences,1,'ascend');
		key = sortrows(key,1,'ascend');
		summ = [key,occurences(:,2)]; % summary per clusterting for occurence
		summ = summ(:,[2 1 3]);
		% add clusterid to overall
		hodor =[zeros(size(table2array(cell2table(udnclustNew)),1),1),table2array(cell2table(udnclustNew))];
		for row = 1:size(hodor,1)
			hodor(row,1)=key(find(strcmp(key(:,1),hodor(row,2))),2);
		end
		hodor = hodor(:,[3 2 1]);
		hodorn = cell(size(hodor,1),3);
		hodorn(:,1)=table2cell(array2table(str2double(hodor(:,1))));
		hodorn(:,2)=udnclustNew(:,1);
		hodorn(:,3)=table2cell(array2table(str2double(hodor(:,3))));
		%%% save on struc 
		eval(strcat('UDN.udnCluster',num2str(cl),'.summary = summ;'));
		eval(strcat('UDN.udnCluster',num2str(cl),'.IdUdnClnmbr = hodorn;'));
		
	end		
	%%% AGAIN FOR ALL CLUSTERS
	x = [clust.cluster1;clust.cluster2;clust.cluster3;clust.cluster4;clust.cluster5;clust.cluster6;clust.cluster7];
	x = x(:,[2 3 4 5 7]);
	%%% IF WANTED WITH TRANSPORT REACTIONS
	%x = [num2cell(deviationMatrix),rxnpathway(:,2)];	
	%%%
	udnclust = cell(size(x,1),3);
	for row = 1:size(x,1)			
		udnclust{row,1} = x{row,2}-x{row,1};
		udnclust{row,2} = x{row,3}-x{row,2};
		udnclust{row,3} = x{row,4}-x{row,3};
	end
	% make changes that are basically 0 to 0
	udnclust = cell2mat(udnclust);
	udnclust(abs(udnclust)<1e-4) = 0;
	udnclust = num2cell(udnclust);
	% exchange with u d n for how the change behaves 
	for row = 1:size(udnclust,1)
		for col = 1:size(udnclust,2)				
			% check for posneg0value
			if udnclust{row,col}>0
				udnclust{row,col} = "U";
			elseif udnclust{row,col}<0
				udnclust{row,col} = "D";
			elseif isequal(udnclust{row,col},0)
				udnclust{row,col} = "N";
			end			
		end
	end		
	% add id
	udnclust = [udnclust,x(:,end)];
	% save udn of that specific cluster 
	UDN.udnClusterAll.udnclust=udnclust;
	%%% tighten frame to classes 
	udnclustNew = cell(size(udnclust,1),2);
	for row = 1:size(udnclust,1)
		udnclustNew{row,1} = strcat(udnclust{row,1},udnclust{row,2},udnclust{row,3});
	end
	udnclustNew(:,2)=udnclust(:,4);
	% save
	UDN.udnClusterAll.udnclust=udnclustNew;	
	%%%
	%%% evaluate udn by 
	[reg,~,ic] = unique(table2array(cell2table(udnclustNew(:,1))));
	occurences = [reg,accumarray(ic,1)];
	gesamt = [	table2array(cell2table(udnclustNew(:,1))),ic];
	%gesamt = [ic,table2array(cell2table(udnclustNew))];
	key = unique(gesamt,'rows');
	% sort
	occurences = sortrows(occurences,1,'ascend');
	key = sortrows(key,1,'ascend');
	summ = [key,occurences(:,2)]; % summary per clusterting for occurence
	summ = summ(:,[2 1 3]);
	% add clusterid to overall
	hodor =[zeros(size(table2array(cell2table(udnclustNew)),1),1),table2array(cell2table(udnclustNew))];
	for row = 1:size(hodor,1)
		hodor(row,1)=key(find(strcmp(key(:,1),hodor(row,2))),2);
	end
	hodor = hodor(:,[3 2 1]);
	hodorn = cell(size(hodor,1),3);
	hodorn(:,1)=table2cell(array2table(str2double(hodor(:,1))));
	hodorn(:,2)=udnclustNew(:,1);
	hodorn(:,3)=table2cell(array2table(str2double(hodor(:,3))));
	%%% 
	UDN.udnClusterAll.summary = summ;
	UDN.udnClusterAll.IdUdnClnmbr = hodorn;
	
	
	
	
	
	
	
	
	
	
	
	
	
end