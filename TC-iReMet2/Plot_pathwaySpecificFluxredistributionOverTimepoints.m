function [] = Plot_pathwaySpecificFluxredistributionOverTimepoints()
    
    %%% read in tables
    % read in rxnlist with corresponding pathways 
    rxnpathway = readtable('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Data/RxnsList_ID_shortname_pathway.csv');
    % read in deviationmatrix of fluxes between WTMT
    deviationMatrix = dlmread('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/day/without_constCofactorRatios/with_slacks/with_slack_minimization/1e8/rep1/deviationMatrix_-_norm_to_Col0.dat','\t');
    maxi = max(max(deviationMatrix));
    
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
        eval(strcat ('pway.p' ,num2str(pathway), ' = [deviationMatrix(find(foundpathway(:,' ,num2str(pathway), ')),:), find(foundpathway(:,' ,num2str(pathway), '))];'));
    end
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%
    %%% PLOTTING %%%
    %%%%%%%%%%%%%%%%
    
    %%% create vertical barplots of flux deviations
    for pathway = 1:size(foundpathway,2)
        eval(strcat('x = pway.p',num2str(pathway),';'));
        bar(x(:,1:4));
        title(strcat('Reaction pathway: ',' ',pathways{pathway}));
        set(gca,'XTickLabel',(rxnpathway(x(:,5),3)'));
        xtickangle(45);
        %ylim([-1 1]);
        %set(gca,'YScale','log');
        ylabel('Distance in [ mol^2 / (gDW d)^2]');
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/pathway_bar_vertical/pathway',num2str(pathway),'_',regexprep(pathways{pathway},'/',''),'.png'));
    end
    
    %%% create horizontal barplots of flux deviations
    for pathway = 1:size(foundpathway,2)
        eval(strcat('x = pway.p',num2str(pathway),';'));
        barh(x(:,1:4));
        %bar(x(:,1:4)/max(max(x)));
        title(strcat('Reaction pathway: ',' ',pathways{pathway}));
        set(gca,'YTickLabel',(rxnpathway(x(:,5),3)'));
        xlabel('Distance in [ mol^2 / (gDW d)^2]');
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/pathway_bar_horizontal/pathway',num2str(pathway),'_',regexprep(pathways{pathway},'/',''),'.png'));
    end
    
    %%% create vertical barplots normalized to maxvalue of flux deviations
    for pathway = 1:size(foundpathway,2)
        eval(strcat('x = pway.p',num2str(pathway),';'));
        bar(x(:,1:4)/max(max(x)));
        title(strcat('Reaction pathway: ',' ',pathways{pathway}));
        set(gca,'XTickLabel',(rxnpathway(x(:,5),3)'));
        xtickangle(45);
        ylabel('Distance in [ mol^2 / (gDW d)^2]');
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/pathway_bar_vertical_maxnormalized/pathway',num2str(pathway),'_',regexprep(pathways{pathway},'/',''),'.png'));
    end
    
    %%% create horizontal barplots normalized to maxvalue of flux deviations
    for pathway = 1:size(foundpathway,2)
        eval(strcat('x = pway.p',num2str(pathway),';'));
        barh(x(:,1:4)/max(max(x)));
        title(strcat('Reaction pathway: ',' ',pathways{pathway}));
        set(gca,'YTickLabel',(rxnpathway(x(:,5),3)'));
        xlabel('Distance in [ mol^2 / (gDW d)^2]');
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/pathway_bar_horizontal_maxnormalized/pathway',num2str(pathway),'_',regexprep(pathways{pathway},'/',''),'.png'));
    end
    
    %{
    %%% heatmap maybe 
    heatmap(deviationMatrix(1:50,1:2));
    normdev = deviationMatrix/max(max(deviationMatrix));
    %normdev = (deviationMatrix-min(min(deviationMatrix)))/(max(max(deviationMatrix))-min(min(deviationMatrix)));
    %normdev = deviationMatrix/(max(max(deviationMatrix))-min(min(deviationMatrix)));
    heatmap(normdev, [], [], '%0.2f', 'Colormap', 'money', ...
        'Colorbar', true);
    h = heatmap(normdev(1:50,1:3));
    h = heatmap(deviationMatrix(1:50,1:3));
    h.Title = 'this is my heatmap lul';
    h.XLabel = 'dayz';
    h.YLabel = 'rxns';
    h.Colormap = 'cool';
    h.ColorMethod = 'mean';
    h.ColorLimits = [-1 1];
    h.GridVisible = 'off';
    %}
    
    %%% plot overall deviation
    overalldev = [5e-02, 1.081e+01, 4.586e+01, 9.086e+01];
    bar((overalldev));
    %title('Distance over all timepoints: 0d-1d-7d-21d');
    set(gca,'XTickLabel',{'day_0','day_1','day_7','day_2_1'},'FontSize',27);
    set(gca,'YScale','log')
	ylim([1e-2 1e3])
    ylabel('Distance [ mol^2 / (gDW d)^2]','FontSize',30);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
	text(1:length(overalldev),overalldev,num2str(overalldev'),'vert','bottom','horiz','center','Fontsize',40); 
    saveas(gcf,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/Overall_Distance.png'));

    %%% plot biom distance 
    bar(deviationMatrix(end,:));
    title('Distance of BiomassReaction over all timepoints: 0d-1d-7d-21d');
    set(gca,'XTickLabel',{'day_0','day_1','day_7','day_2_1'});
    ylabel('Distance in [ mol^2 / (gDW d)^2]','FontSize',40);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/Overall_Distance_biomass.png'));

    
    %%% plot fluxbehaviour over timepoints as lines 
    for pathway = 1:size(foundpathway,2)
        
        % close and open new figure 
        close(gcf);
        figure();
        
        % initialize important plot parameters
        xticks([1:4]);
        curvespecs = {'r-','g-','b-','y-','c-','r--','g--','b--','y--','c--',...
                    'r-.','g-.','b-.','y-.','c-.','r:','g:','b:','y:','c:',};        
        % assign pathway matrix
        eval(strcat('x = pway.p',num2str(pathway),';'));
        disp(size(x));
        % break if it is too large
        if size(x,1)>length(curvespecs)
            break
        end        
                
        % open plot and set axis
        hold on;
        set(gca,'xticklabel',{'day_0','day_1','day_7','day_2_1'});
        
        % vector for curvenames
        legi=[];
        % plot all curves 
        for i = 1:size(x,1)
            plot(x(i,1:4),curvespecs{i},'LineWidth',1);
            legi = [legi,''',''', rxnpathway{x(i,5),3}];
        end
        legi(1:3) = [];        
        eval(strcat('names = {''',legi,''',''0 Line''};'));                
        % add horizontal line at 0 
        y = [0,0,0,0];
        plot(y,'k-','LineWidth',2);
        
        %plot legend etc 
        legend(names,'Location','northwestoutside'); 
        title(strcat('Reaction pathway: ',' ',pathways{pathway}));
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/pathway_linebehaviour/pathway',num2str(pathway),'_',regexprep(rxnpathway{x(1,5),4},'/',''),'.png'));
        
        hold off;        
    end
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% create vertical barplots of flux deviations
    for pathway = 1:size(foundpathway,2)
        eval(strcat('x = pway.p',num2str(pathway),';'));
        bar(x(:,1:4)./maxi(1:4));
        title(strcat('Reaction pathway: ',' ',pathways{pathway}));
        set(gca,'XTickLabel',(rxnpathway(x(:,5),3)'));
        xtickangle(45);
        ylim([-1 1]);
        ylabel('Distance in [ mol^2 / (gDW d)^2]');
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf,strcat('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/pathway_bar_vertical_normalizedToOverallMax/pathway',num2str(pathway),'_',regexprep(pathways{pathway},'/',''),'.png'));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%
    %%% PLAYGROUND %%%
    %%%%%%%%%%%%%%%%%%
    
    data = deviationMatrix;
    WT_fluxmatrix = dlmread('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/day/without_constCofactorRatios/with_slacks/with_slack_minimization/1e8/rep1/wildtypeFluxMatrix_-_norm_to_Col0.dat','\t');
    MT_fluxmatrix = dlmread('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/day/without_constCofactorRatios/with_slacks/with_slack_minimization/1e8/rep1/mutantFluxMatrix_-_norm_to_Col0.dat','\t');
    WT = WT_fluxmatrix;
    MT = MT_fluxmatrix;
    
    
    % corr matrix
    mpear = zeros(4);
    for i = 1:4
        for j = 1:4
            mpear(i,j)=corr(data(:,i),data(:,j));
        end
    end
    
    % welche reaction ist die von tag 1 auf tag 2 ab staerksten
    % unterschiedliche 
    gradient = [data(:,1),data(:,2)];
    t1t2 = sqrt(((data(:,1)-data(:,2)).^2));
    [t1t2s,t] = sort(t1t2,'descend');
    rxnpathway(t(1:20),:)
    
    
    
    %%% read in idex table 
    idrxn = table2array(readtable('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/toprxnIDs.csv'));
    idrxn(:,3)
    
    % 0d topx% reactions fluxes and names
    top0 = [num2cell(data(idrxn(:,3),1)),rxnpathway(idrxn(:,3),2),rxnpathway(idrxn(:,3),3),rxnpathway(idrxn(:,3),4)];
    % 1d topx% reactions fluxes and names
    top1 = [num2cell(data(idrxn(:,6),1)),rxnpathway(idrxn(:,6),2),rxnpathway(idrxn(:,6),3),rxnpathway(idrxn(:,6),4)];
    % 7d topx% reactions fluxes and names
    top7 = [num2cell(data(idrxn(:,9),1)),rxnpathway(idrxn(:,9),2),rxnpathway(idrxn(:,9),3),rxnpathway(idrxn(:,9),4)];
    % 21d topx% reactions fluxes and names
    top21 = [num2cell(data(idrxn(:,12),1)),rxnpathway(idrxn(:,12),2),rxnpathway(idrxn(:,12),3),rxnpathway(idrxn(:,12),4)];
    
    [d,t]=sort(data(:,3),'descend')
    rxnpathway(t,:)
    
    
    %%% create table of overall top regulated reactions per timepoint
    d1 = data(:,1);
    [d1s,id]=sort(d1,'descend');
    t0 = [num2cell(d1s),rxnpathway(id,2:4),num2cell(id)];
    %
    d2 = data(:,2);
    [d2s,id]=sort(d2,'descend');
    t1 = [num2cell(d2s),rxnpathway(id,2:4),num2cell(id)];
    %
    d3 = data(:,3);
    [d3s,id]=sort(d3,'descend');
    t2 = [num2cell(d3s),rxnpathway(id,2:4),num2cell(id)];
    %
    d4 = data(:,4);
    [d4s,id]=sort(d4,'descend');
    t3 = [num2cell(d4s),rxnpathway(id,2:4),num2cell(id)];
    
    name = {'fluxdistance','abbreviation','rxnName','pathway','rxnID'};
    writetable(cell2table(t0,'VariableNames',name),'/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/allrxnsorted_t0.csv');
    writetable(cell2table(t1,'VariableNames',name),'/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/allrxnsorted_t1.csv');
    writetable(cell2table(t2,'VariableNames',name),'/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/allrxnsorted_t2.csv');
    writetable(cell2table(t3,'VariableNames',name),'/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/figures/allrxnsorted_t3.csv');
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % corr wts to mts
    % create giant table 
    C = [WT,MT];
    
    % show me the reactions of WT or MT that are different from 0d to 1d
    rxnpathway(find((WT(:,2)-WT(:,1))>0),:)
    rxnpathway(find((MT(:,2)-MT(:,1))>0),:)
    
    % show me the reactions of WT or MT that are different from 1d to 7d
    rxnpathway(find((WT(:,3)-WT(:,2))>0),:)
    rxnpathway(find((MT(:,3)-MT(:,2))>0),:)
    
    % show me the reactions of WT or MT that are different from 7d to 21d
    rxnpathway(find((WT(:,4)-WT(:,3))>0),:)
    rxnpathway(find((MT(:,4)-MT(:,3))>0),:)
 
    %%% euclidian distance of type at t_n to t_n+1
    disWT0dto1d = sqrt(sum((WT(:,1) - WT(:,2)) .^ 2));
    disWT1dto7d = sqrt(sum((WT(:,2) - WT(:,3)) .^ 2));
    disWT7dto21d = sqrt(sum((WT(:,3) - WT(:,4)) .^ 2));
    disMT0dto1d = sqrt(sum((MT(:,1) - MT(:,2)) .^ 2));
    disMT1dto7d = sqrt(sum((MT(:,2) - MT(:,3)) .^ 2));
    disMT7dto21d = sqrt(sum((MT(:,3) - MT(:,4)) .^ 2));
    
    disWT0dtoMT1d = sqrt(sum((WT(:,1) - MT(:,2)) .^ 2));
    disWT1dtoMT7d = sqrt(sum((WT(:,2) - MT(:,3)) .^ 2));
    disWT7dtoMT21d = sqrt(sum((WT(:,3) - MT(:,4)) .^ 2));
    
    disWT1dtoMT0d = sqrt(sum((WT(:,2) - MT(:,1)) .^ 2));
    disWT7dtoMT1d = sqrt(sum((WT(:,3) - MT(:,2)) .^ 2));
    disWT21dtoMT7d = sqrt(sum((WT(:,4) - MT(:,3)) .^ 2));
    
    a = [ disWT0dto1d,disWT1dto7d,disWT7dto21d];
    b = [ disMT0dto1d,disMT1dto7d,disMT7dto21d];
    c = [ disWT0dtoMT1d,disWT1dtoMT7d,disWT7dtoMT21d];
    d = [ disWT1dtoMT0d,disWT7dtoMT1d,disWT21dtoMT7d];
    n = [ a;b;c;d];
    
    %%% 
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% create heatmap test
	x = pway.p2(:,1:4);
	xb = pway.p2;
	n = num2cell(xb);
	%%% exchange with name %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	%
	% CORRECT NAMING PART OF THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	%
	for i=1:size(n,1)		
		id = n{i,5};
		n{i,5} = rxnpathway{i,2};		
	end 	
	xt = cell2table(n,'VariableNames',{'t0','t1','t2','t3','nmbr'});
	%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	%%% transform xt
	% empty thing 
	c = 1;
	
	%rxnnamesc = ones(1,length(x(:)))*NaN;
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
	h.FontSize = 12;
	%h.ColorLimits = [-0.003 0.003];
	%h.ColorLimits = [-1.5 1.5];
	h.ColorLimits = [-5 5];
	%h.ColorLimits = [-10 10];
	colormap(redblue())

	
	hc = clustergram(x);
	
	
	
	
	
    % make just values show 
    fig = figure;
	colormap('hot')
	imagesc(x)
	
	%%%%%%%%%%%%%%%%
    %%% HEATMAP GESAMTBILD THINGY
	%%%%%%%%%%%%%%%%
	
	data = deviationMatrix;
	data(abs(data)<1e-4) = 0;
	d = data;
	d = [d,(1:549)'];
	x = d;
	
	
	
	
	
    


end 