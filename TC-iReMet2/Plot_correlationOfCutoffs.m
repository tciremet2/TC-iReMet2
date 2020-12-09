function [] = Plot_correlationOfCutoffs(rvCorrStruct,AnalysisParameters,days,maxratio_array)

    %disp(rvCorrStruct);
    %disp(kNumbRxns);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CALCULATE CORRELATION ARRAYS %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % init rvCorrStruct
    for cut=1:length(maxratio_array)
        for time=0:days-1
            eval(strcat('rvCorrStruct.e',num2str(cut),'.t',num2str(time),' = [rvCorrStruct.e',num2str(cut),'.main_WT.t',num2str(time),'; rvCorrStruct.e',num2str(cut),'.main_MT.t',num2str(time),'];' ));
            %disp(strcat('rvCorrStruct.e',num2str(cut),'.t',num2str(time),' = [rvCorrStruct.e',num2str(cut),'.main_WT.t',num2str(time),'; rvCorrStruct.e',num2str(cut),'.main_MT.t',num2str(time),']' ));
        end
    end
    
    % init correlation vectors
    for i=0:days-1
        eval(strcat('t',num2str(i),'=zeros(1,length(maxratio_array)-1);'));
        %disp(strcat('t',num2str(i),'=zeros(1,length(maxratio_array)-1);'))
    end
    
    % calculate correlations as rvcoeff
    for day=0:days-1
        for cor=1:length(maxratio_array)-1
            eval(strcat('t',num2str(day),'(',num2str(cor),') = corr2( rvCorrStruct.e',num2str(cor),'.t',num2str(day),',rvCorrStruct.e',num2str(cor+1),'.t',num2str(day),');'))
            %disp(strcat('t',num2str(day),'(',num2str(cor),') = corr2( rvCorrStruct.e',num2str(cor),'.t',num2str(day),',rvCorrStruct.e',num2str(cor+1),'.t',num2str(day),');'));
        end
    end
    % calculate correlations as vector with just pearson 
    for day=0:days-1
        for cor=1:length(maxratio_array)-1
            eval(strcat('t',num2str(day),'(',num2str(cor),') = corr( rvCorrStruct.e',num2str(cor),'.t',num2str(day),'(:,1),rvCorrStruct.e',num2str(cor+1),'.t',num2str(day),'(:,1));'));
            %disp(strcat('t',num2str(day),'(',num2str(cor),') = corr( rvCorrStruct.e',num2str(cor),'.t',num2str(day),'(:,1),rvCorrStruct.e',num2str(cor+1),'.t',num2str(day),'(:,1));'));
        end
    end

    % for all array thingy
    tall = zeros(1,length(maxratio_array)-1);
    for cor=1:length(maxratio_array)
        % add to one table - hardcoded change this later 
        %disp(strcat('rvCorrStruct.all.e',num2str(cor),'= [rvCorrStruct.e',num2str(cor),'.t0',';rvCorrStruct.e',num2str(cor),'.t1',';rvCorrStruct.e',num2str(cor),'.t2',';rvCorrStruct.e',num2str(cor),'.t3];' ));
        eval(strcat('rvCorrStruct.all.e',num2str(cor),'= [rvCorrStruct.e',num2str(cor),'.t0',';rvCorrStruct.e',num2str(cor),'.t1',';rvCorrStruct.e',num2str(cor),'.t2',';rvCorrStruct.e',num2str(cor),'.t3];' ));
    end
    
    for cor=1:length(maxratio_array)-1
        %disp(strcat( 'tall(',num2str(cor),')= corr2( rvCorrStruct.all.e' ,num2str(cor),',rvCorrStruct.all.e',num2str(cor+1),');'));
        eval(strcat( 'tall(',num2str(cor),')= corr2( rvCorrStruct.all.e' ,num2str(cor),',rvCorrStruct.all.e',num2str(cor+1),');'));
    end   
    
    % for all array thingy for vector version 
    for cor=1:length(maxratio_array)-1
        %disp(strcat( 'tall(',num2str(cor),')= corr( rvCorrStruct.all.e' ,num2str(cor),'(:,1),rvCorrStruct.all.e',num2str(cor+1),'(:,1));'));
        eval(strcat( 'tall(',num2str(cor),')= corr( rvCorrStruct.all.e' ,num2str(cor),'(:,1),rvCorrStruct.all.e',num2str(cor+1),'(:,1));'));
    end  


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%
    %%% PLOT %%% 
    %%%%%%%%%%%%
    
    % assign time specific vectors
    a = t0;     % time point 0 -- 0d
    b = t1;     % time point 1 -- 1d
    c = t2;     % time point 2 -- 7d
    d = t3;     % time point 3 -- 21d
    e = tall;   % time point all_timepoints -- all_d
  
    name = {};
    for i=1:length(maxratio_array)
        eval(strcat('name{',num2str(i),'}="e',num2str(i),'-',num2str(i+1),'";'));
        %disp(strcat('name{',num2str(i),'}="e',num2str(i),'"'));
    end

    % build figure
    figure('Renderer', 'painters', 'Position', [10 10 1000 500])
    xlabel('maxratio as e^x')
    ylabel('rv-coeff')
    title('correlation of fluxes from e^1 to e^n^+^1')
    hold on;
    plot(a,'b-','LineWidth',2)
    xticks([1:21])
    xticklabels = name;
    set(gca,'xticklabel',name)
    plot(b,'r-','LineWidth',2)
    plot(c,'g-','LineWidth',2)
    plot(d,'m-','LineWidth',2)
    plot(e,'c--','LineWidth',2)

    %legend({'a = 0 day','b = 1 day','c = 7 day', 'd = 21 day'},'Location','northwest')
    legend({'a = day 0','b = day 1','c = day 7', 'd = day 21', 'e = all day'})

    %set(gca,'FontSize',8)
    
    % save figure
    saveas(gcf,'/home/mpimp-golm.mpg.de/pries1347/winhome/main/MAIN-uebergabe/TC-iReMet2/Results/figures/CorrCutOffs.jpg')

    
    
end