function [summDeviationMatrix] = Eval_enrichment()

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% READ IN TABLES AND CREATE PATHWAY ARRAYS/INFORMATION %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%% clear workspace
	clear;
	%%% read in data and set important variables
	kNumbRxns = 549;
	rxnpathway = readtable('E:\Uni\MasterarbeitBzwPaper\TCiReMet2_all\MASTERARBEIT\MAIN-uebergabe\TC-iReMet2\Data\RxnsList_ID_shortname_pathway.csv');
    % read in deviationmatrix of fluxes between WTMT
    deviationMatrix = dlmread('E:\Uni\MasterarbeitBzwPaper\TCiReMet2_all\MASTERARBEIT\MAIN-uebergabe\TC-iReMet2\Results\day\without_constCofactorRatios\with_slacks\with_slack_minimization\1e8\rep1\deviationMatrix_-_norm_to_Col0.dat','\t');
    deviationMatrix(abs(deviationMatrix)<1e-4) = 0;
	load artificalRxns.mat
	% convert to cell for easier handling 
    rxnpathway = table2cell(rxnpathway);
    % extract pathways and exclude empty cells 
    pathways = rxnpathway(:,5);
    pathways = pathways(~cellfun('isempty',pathways));
	pathways = unique(pathways);
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% CREATE PATHWAY ARRAYS WITH TRANSPORT REACTIONS %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%% create array where 1 corresponds to present reaction in pathway(column) 
    foundpathway = ones(size(rxnpathway,1),size(pathways,1))*NaN;
    for pathway = 1:size(pathways,1)
        foundpathway(:,pathway)=contains(rxnpathway(:,4),pathways{pathway});
	end

	%%% create table of fluxdeviation and reactionIDs for plotting 
    for pathway = 1:size(foundpathway,2)
        eval(strcat ('pway.withTnA.p' ,num2str(pathway), '.M = [deviationMatrix(find(foundpathway(:,' ,num2str(pathway), ')),:), find(foundpathway(:,' ,num2str(pathway), '))];'));
        eval(strcat ('pway.withTnA.p' ,num2str(pathway), '.pwayname = pathways(', num2str(pathway), ');'));
	end
	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% CREATE PATHWAY ARRAYS WITHOUT TRANSPORT REACTIONS %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%% exclude artifical and transport reactions 
	foundpathway = foundpathway(setdiff(1:549,artificial_rxns_and_transporters)',:);
	rxnpathway = rxnpathway(setdiff(1:549,artificial_rxns_and_transporters)',:);
    
	%%% create array where 1 corresponds to present reaction in pathway(column)
	for pathway = 1:size(pathways,1)
        foundpathway(:,pathway)=contains(rxnpathway(:,4),pathways{pathway});
	end
	%%% create table of fluxdeviation and reactionIDs for plotting 
    for pathway = 1:size(foundpathway,2)
        eval(strcat ('pway.withoutTnA.p' ,num2str(pathway), '.M = [deviationMatrix(find(foundpathway(:,' ,num2str(pathway), ')),:), find(foundpathway(:,' ,num2str(pathway), '))];'));
        eval(strcat ('pway.withoutTnA.p' ,num2str(pathway), '.pwayname = pathways(', num2str(pathway), ');'));
	end

	%%% revert to original
	rxnpathway = readtable('E:\Uni\MasterarbeitBzwPaper\TCiReMet2_all\MASTERARBEIT\MAIN-uebergabe\TC-iReMet2\Data\RxnsList_ID_shortname_pathway.csv');
	rxnpathway = table2cell(rxnpathway);
    pathways = rxnpathway(:,5);
    pathways = pathways(~cellfun('isempty',pathways));
	pathways = unique(pathways);
    foundpathway = ones(size(rxnpathway,1),size(pathways,1))*NaN;
    for pathway = 1:size(pathways,1)
        foundpathway(:,pathway)=contains(rxnpathway(:,4),pathways{pathway});
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% DESCRIPTIVE STATISTICS OF FLUX DEVIATIONS %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	summDeviationMatrix = summary(array2table(deviationMatrix));
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% RANK REACTIONS FOR ALL TIME POINTS AND BETWEEN EVERY TIME POINT STEP %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%% create initial based tables 
	nominalDeviationmatrix = deviationMatrix;
	relativeDeviationmatrix = deviationMatrix./max(abs(deviationMatrix),[],2);
	relativeDeviationmatrix(isnan(relativeDeviationmatrix)) = 0;
	nominalRegulationvalueTimepoint = ones(size(deviationMatrix,1), size(deviationMatrix,2))*NaN; % col1=t1-t2,col2=t2-t3,col3=t3-t4,col4=overall
	relativeRegulationvalueTimepoint = ones(size(deviationMatrix,1), size(deviationMatrix,2))*NaN; % col1=t1-t2,col2=t2-t3,col3=t3-t4,col4=overall
	
	%%% calculate regulation values for nominal and relative
	nominalRegulationvalueTimepoint(:,1) = abs(nominalDeviationmatrix(:,1)-nominalDeviationmatrix(:,2));
	nominalRegulationvalueTimepoint(:,2) = abs(nominalDeviationmatrix(:,2)-nominalDeviationmatrix(:,3));
	nominalRegulationvalueTimepoint(:,3) = abs(nominalDeviationmatrix(:,3)-nominalDeviationmatrix(:,4));
	nominalRegulationvalueTimepoint(:,4) = nominalRegulationvalueTimepoint(:,1)+nominalRegulationvalueTimepoint(:,2)+nominalRegulationvalueTimepoint(:,3);	
	%
	relativeRegulationvalueTimepoint(:,1) = abs(relativeDeviationmatrix(:,1)-relativeDeviationmatrix(:,2));
	relativeRegulationvalueTimepoint(:,2) = abs(relativeDeviationmatrix(:,2)-relativeDeviationmatrix(:,3));
	relativeRegulationvalueTimepoint(:,3) = abs(relativeDeviationmatrix(:,3)-relativeDeviationmatrix(:,4));
	relativeRegulationvalueTimepoint(:,4) = relativeRegulationvalueTimepoint(:,1)+relativeRegulationvalueTimepoint(:,2)+relativeRegulationvalueTimepoint(:,3);	
	%
	nominalRegulationvalueTimepoint = [nominalRegulationvalueTimepoint,[1:549]'];
	relativeRegulationvalueTimepoint = [relativeRegulationvalueTimepoint,[1:549]'];
	nominalRegulationvalueTimepoint = nominalRegulationvalueTimepoint(setdiff(1:549,artificial_rxns_and_transporters)',:);
	relativeRegulationvalueTimepoint = relativeRegulationvalueTimepoint(setdiff(1:549,artificial_rxns_and_transporters)',:);
	%
	summNom = summary(array2table(nominalRegulationvalueTimepoint));
	summRel = summary(array2table(relativeRegulationvalueTimepoint));
	
	
	%%% init cell of significance testing results for both measures
	nominalSig = cell(size(pathways,1),4);
	relativeSig = cell(size(pathways,1),4);
	
	%%% sets of rxn ids for reg or not at each time step
	regulation = [];
	for time = 1:size(deviationMatrix,2)
		
		% set median as cutoffs
		eval(strcat('medianNom = summNom.nominalRegulationvalueTimepoint' ,num2str(time), '.Median;'));
		eval(strcat('medianRel = summRel.relativeRegulationvalueTimepoint' ,num2str(time), '.Median;')); 
		
		%%% make sets of rxn ids 
		regulatedNom = nominalRegulationvalueTimepoint(nominalRegulationvalueTimepoint(:,time)>medianNom,5);
		nonregulatedNom = nominalRegulationvalueTimepoint(nominalRegulationvalueTimepoint(:,time)<=medianNom,5);
		%	
		regulatedRel = relativeRegulationvalueTimepoint(relativeRegulationvalueTimepoint(:,time)>medianRel,5);
		nonregulatedRel = relativeRegulationvalueTimepoint(relativeRegulationvalueTimepoint(:,time)<=medianRel,5);
		
		%%% save on struc 
		eval(strcat('regulation.nominal.x', num2str(time), '.regulated = regulatedNom;'));
		eval(strcat('regulation.nominal.x', num2str(time), '.nonregulated = nonregulatedNom;'));
		%
		eval(strcat('regulation.relative.x', num2str(time), '.regulated = regulatedRel;'));
		eval(strcat('regulation.relative.x', num2str(time), '.nonregulated = nonregulatedRel;'));
		
	end
	
	%%% add info about pvalue and odds ratio
	nominalPvalue = cell(size(pathways,1),4);
	relativePvalue = cell(size(pathways,1),4);
	nominalOdds = cell(size(pathways,1),4);
	relativeOdds = cell(size(pathways,1),4);
	
	%%%
	for type = 1:2
		if type == 1
			s = 'nominal';
		elseif type == 2 
			s = 'relative';
		end
		
		for timepoint = 1:size(deviationMatrix,2)
			
			eval(strcat('regulatedrxns = regulation.', s, '.x', num2str(timepoint), '.regulated;'));
			kNumbregulatedrxns = size(regulatedrxns,1);
			eval(strcat('nonregulatedrxns = regulation.', s, '.x', num2str(timepoint) ,'.nonregulated;'));
			kNumbnonregulatedrxns = size(regulatedrxns,1);
			
			for pathway = 1:size(pathways,1)
				
				eval(strcat('x = pway.withoutTnA.p',num2str(pathway),'.M;'));
				
				numbRegRxnsInPathway = length(intersect(x(:,end),regulatedrxns));
				numbnonRegRxnsInPathway = size(x,1)-numbRegRxnsInPathway;
				
				%contingencyTable = [numbRegRxnsInPathway,numbnonRegRxnsInPathway;
				%					(kNumbregulatedrxns-numbRegRxnsInPathway),(kNumbnonregulatedrxns-numbnonRegRxnsInPathway)];
				
				contingencyTable = [(kNumbnonregulatedrxns-numbnonRegRxnsInPathway),numbnonRegRxnsInPathway;
									(kNumbregulatedrxns-numbRegRxnsInPathway),numbRegRxnsInPathway];
								
				[fish,p,stats] = fishertest(contingencyTable,'Tail','right','Alpha',0.01); % maybe retrieve actual p value here 
				
				eval(strcat(s, 'Sig{pathway,timepoint} = fish;'));
				eval(strcat(s, 'Pvalue{pathway,timepoint} = p;'));
				eval(strcat(s, 'Odds{pathway,timepoint} = stats.OddsRatio;'));
			end
		end
	end		
	
	%%% save ablauf 
	pathways(table2array(cell2table(relativeSig(:,4))));
	find(table2array(cell2table(relativeSig(:,4))));
	relativePvalue(find(table2array(cell2table(relativeSig(:,4)))),4);
	relativeOdds(find(table2array(cell2table(relativeSig(:,4)))),4);
	
	% dont do this before below if below needed 
	relativeSig2 = [cell2table(pathways),cell2table(relativeSig)];
	nominalSig2 = [cell2table(pathways),cell2table(nominalSig)];
	writetable(relativeSig2,'//home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MAIN-uebergabe/TC-iReMet2/Results/relativeSig.csv');
	writetable(nominalSig2,'/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MAIN-uebergabe/TC-iReMet2/Results/nominalSig.csv');
	%	
	allrel=[pathways(table2array(cell2table(relativeSig(:,4)))),relativePvalue(find(table2array(cell2table(relativeSig(:,4)))),4),relativeOdds(find(table2array(cell2table(relativeSig(:,4)))),4)];
	allnom=[pathways(table2array(cell2table(nominalSig(:,4)))),nominalPvalue(find(table2array(cell2table(nominalSig(:,4)))),4),nominalOdds(find(table2array(cell2table(nominalSig(:,4)))),4)];
	% for 0dto1d
	allrel01=[pathways(table2array(cell2table(relativeSig(:,1)))),relativePvalue(find(table2array(cell2table(relativeSig(:,1)))),1),relativeOdds(find(table2array(cell2table(relativeSig(:,1)))),1)];
	allnom01=[pathways(table2array(cell2table(nominalSig(:,1)))),nominalPvalue(find(table2array(cell2table(nominalSig(:,1)))),1),nominalOdds(find(table2array(cell2table(nominalSig(:,1)))),1)];
	% save stuff
	writetable(cell2table(allrel),strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MAIN-uebergabe/TC-iReMet2/Results/allrel.csv'));
	writetable(cell2table(allnom),strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MAIN-uebergabe/TC-iReMet2/Results/allnom.csv'));
	writetable(cell2table(allrel01),strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MAIN-uebergabe/TC-iReMet2/Results/allrel01.csv'));
	writetable(cell2table(allnom01),strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MAIN-uebergabe/TC-iReMet2/Results/allnom01.csv'));

	% early timepoints
	rxnpathway(regulation.relative.x1.regulated(1:3,:),1:4);
	rxnpathway(regulation.nominal.x1.regulated(1:3,:),1:4);
	
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% ENRICHMENT CLUSTERS %%% 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%{
	% read in cluster data
	clust = [];
	for i = 1:7	
		eval(strcat('clust.cluster',num2str(i),' = readtable(''/home/mpimp-golm.mpg.de/pries1347/winhome/main/MAIN-uebergabe/TC-iReMet2/Results/R_calcs/cluster', num2str(i), '.txt'');'));
		eval(strcat('clust.cluster',num2str(i),'=table2cell(clust.cluster',num2str(i),');'));
	end
	
	% calculate and updowntable based on size of frame
	UDN = [];
	% evalution 
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
		eval(strcat('UDN.udnCluster',num2str(cl),'.size=size(UDN.udnCluster',num2str(cl),'.IdUdnClnmbr,1);'));
		
		%%% enriched or %tage 
		tu = str2double(summ(:,[1 3]));
		tucount = sum(tu(:,2));
		per = tu(:,2)./tucount;
		eval(strcat('UDN.udnCluster',num2str(cl),'.summary = [summ,string(per)];'));
		
		
		
		
	end	
	
	for cl = 1:7
		eval(strcat('i=max(str2double(UDN.udnCluster',num2str(cl),'.summary(:,end)));'));
		fprintf('%d\n',i)
	end
	
	% find all observed classes, dann fuer jedes cluster wie viel % vom cl
	% drin
	gigaanteile =[];
	for cl = 1:7
		eval(strcat('i=table2array(cell2table(UDN.udnCluster',num2str(cl),'.IdUdnClnmbr(:,2)));'));
		%disp(i)
		gigaanteile = [gigaanteile;i];
	end
	
	
	[reg,~,id] = unique(gigaanteile);
	occurences = [reg,accumarray(id,1)];
	reguclassesfound = size(reg,1);
	
	
	gig = cell(reguclassesfound,8);
	gig(:,1)= table2cell(array2table(occurences(:,1)));
	clsize = [];
	for cl = 1:7
		eval(strcat('i=UDN.udnCluster',num2str(cl),'.summary;'));
		clsize = [clsize,sum(str2double(i(:,3)))];
		
		for row = 1:size(gig,1)
			
			for code = 1:size(i,1)
				
				if gig{row,1} == i(code,2)
					
					gig{row,cl+1} = i(code,3);
					
				end
				
			end
			
		end
				
		
	end
	
	gig(find(cellfun(@isempty,gig))) = {"0"};
	for j = 18:length(gig(:))
		gig{j} = str2double(gig{j});
	end
	
	
	gigadd = cell(size(gig,1),2);
	gig = [gig,gigadd];
	gigadd2 = table2cell(array2table(zeros(1,size(gig,2))));
	gig = [gig;gigadd2];
	%gig = cell2table(gig);
	gig{end,1}="clSize";
	
	for col = 2:8		
		gig{end,col} = sum(table2array(cell2table(gig(:,col))));		
		for row = 1:18			
			gig{row,9} = sum(table2array(cell2table(gig(row,2:8))));			
		end
	end
	for row = 1:17
		gig{row,end}=gig{row,(end-1)}/gig{end,(end-1)};
	end
	
	
	
	
	for cl = 1:7
		eval(strcat('x=UDN.udnCluster',num2str(cl),'.size;'));
		fprintf('%d\n',x)
	end
	%}
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%	make giant significance table %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%% giant for relative change
	SupRelativeSignificance = [pathways,relativePvalue(:,4),relativeOdds(:,4)];
	%%% giant for relative change
	SupNominalSignificance = [pathways,nominalPvalue(:,4),nominalOdds(:,4)];
	%%% save as tables
	writetable(cell2table(SupRelativeSignificance),strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MAIN-uebergabe/TC-iReMet2/Results/SupRelativeSignificance.csv'));
	writetable(cell2table(SupNominalSignificance),strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MAIN-uebergabe/TC-iReMet2/Results/SupNominalSignificance.csv'));
	
	
	
	
	
	
	
	
	
	
	

end