function [] = Plot_AndEvalHeatmaps()

	%%%%%%%%%%%%%%%%%%%%
	%%% READ IN DATA %%%
	%%%%%%%%%%%%%%%%%%%%
	
	% read in rxnlist with corresponding pathways 
	rxnpathway = readtable('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Data/RxnsList_ID_shortname_pathway.csv');
	rxnpathway = table2cell(rxnpathway);
	% read in deviationmatrix of fluxes between WTMT
	deviationMatrix = dlmread('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/day/without_constCofactorRatios/with_slacks/with_slack_minimization/1e8/rep1/deviationMatrix_-_norm_to_Col0.dat','\t');
	% load artificial reactions
	load artificalRxns.mat
	
	%%%%%%%%%%%%%%%%%%%%%%%
	%%% PREPROCESS DATA %%%
	%%%%%%%%%%%%%%%%%%%%%%%
	
	% adjust small values
	data = deviationMatrix;
	data(abs(data)<1e-4) = 0;
	normData = data;
	
	% normalize per row to get behaviour over time for clustering
	for row= 1:size(normData,1)
		% see if row is simply 0
		checker = sum(normData(row,:));
		if ~isequal(checker,0)
			
			normData(row,:) = normData(row,:)./max(abs(normData(row,:)));
		end
	end
		
	% exclude artificial and transport reactions 
	normDataExcl = normData(setdiff(1:549,artificial_rxns_and_transporters)',:);
	% exclude further 0 rows 
	normDataExclN0 = normDataExcl(find(sum(normDataExcl,2)),:);

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% DO KMEANS CLUSTERING %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% init overlord plotstruc 
	clusterH = [];
	% km is the amount of clusters 
	for km = 1:10 
		
		%normed data
		eval(strcat('clusterH.normData.km', num2str(km), '.ID = kmeans(normData,', num2str(km), ');'));
		eval(strcat('clusterH.normData.km', num2str(km), '.M  = [kmeans(normData,', num2str(km), '),normData];'));
		
		% normDataExcl
		eval(strcat('clusterH.normDataExcl.km', num2str(km), '.ID = kmeans(normDataExcl,', num2str(km), ' );'));
		eval(strcat('clusterH.normDataExcl.km', num2str(km), '.M  = [kmeans(normDataExcl,', num2str(km), '),normDataExcl];'));
		
		% normDataExcl
		eval(strcat('clusterH.normDataExclN0.km', num2str(km), '.ID = kmeans(normDataExclN0,', num2str(km), ');'));
		eval(strcat('clusterH.normDataExclN0.km', num2str(km), '.M  = [kmeans(normDataExclN0,', num2str(km), '),normDataExclN0];'));
		
	end
	
	save('test.mat');
	clear
	load test
	
	%%%%%%%%%%%%%%%%%%%%%
	%%% PLOT HEATMAPS %%%
	%%%%%%%%%%%%%%%%%%%%%
	clusterNR = 1;
	x = clusterH.normData.km7.M;
	x = [x,[1:549]'];
	x(:,2:5) = data(:,1:4);
	x = sortrows(x,1,'ascend');
	%
	x = x(x(:,1)==clusterNR,:);
	%
	xb = x(:,2:end);
	n = num2cell(xb);
	for i = 1:size(n,1)
		id = n{i,5};
		n{i,5} = rxnpathway{id,2};
	end
	xt = cell2table(n,'VariableNames',{'t0','t1','t2','t3','nmbr'});
	x = clusterH.normData.km5.M;
	x = [x,[1:549]'];
	x(:,2:5) = data(:,1:4);
	%
	x = x(x(:,1)==clusterNR,:);
	%
	x = x(:,2:5);
	
	% empty thing 
	c = 1;	
	rxnnamesc = strings(1,length(x(:)));
	timepointc = strings(1,length(x(:)));
	valuec = ones(1,length(x(:)))*NaN;
	for row = 1:size(x,1)
		for col = 1:size(x,2)			
			rxnnamesc(c) = n{row,5};
			timepointc(c) = string(xt.Properties.VariableNames{col});
			valuec(c) = x(row,col);
			% update c			
			c = c + 1;		
		end
	end
	% create table out of it 
	T = table(rxnnamesc(:),timepointc(:),valuec(:),'VariableNames',{'reaction','timepoint','fluxdiff'});
	h = heatmap(T,'timepoint','reaction','ColorVariable','fluxdiff');
	h.Title = 'pathway X';
	h.XLabel = 'timepoints';
	h.FontSize = 10;
	%h.ColorLimits = [-0.003 0.003];
	h.ColorLimits = [-10 10];
	%h.ColorLimits = [-5 5];
	%h.ColorLimits = [-10 10];
	colormap(redblue())
	
	
	
	
	
	
	
	
	
	
	
	set(h, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
	saveas(h,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/heatmaps_loop/hcl',  ,'.png'));

		
	
	
	
	
	
	
	
	
	
	
	




end