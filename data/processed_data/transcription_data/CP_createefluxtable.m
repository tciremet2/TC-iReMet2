%%%% CP run eflux on Transcriptomic data table 

%%% read data 
tdata = readtable('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/data_final_tables/transcription_data/final_aggregated_tdata.txt','Delimiter',',');
load('Model_A');
model = Model_A; 
model = Correct_model_mistakes(model);


%%% init
% atgcodelist - list of measured atg codes derived from the read in data
atgcodelist = table2array(tdata(:,1));
% description - stores description of the read in table headers
description = tdata.Properties.VariableNames;
% originaldata - stores the original read in data values 
originaldata = table2array(tdata(:,2:end));
originaldata(:,6)=NaN; % NaN for the empty column 
% eflux_array - stores by eflux calculated values for every experiment  
eflux_array = zeros(size(model.ub,1),size(originaldata,2))*NaN;


%%% calculate eflux values for every experiment/measurement from tdata 
for measurement = 1:size(eflux_array,2)
    eflux_array(:,measurement) = CPeflux_final(model, 1, 1, 0, 0, tdata(:,[1 measurement+1]),1);
end

% get rid of introduced 0s by eflux, there shouldnt be any values possible
% since there are no measurements 
eflux_array(:,6)=NaN;

%%% 
efluxtable = array2table(eflux_array,'Variablenames',description(2:end));


%%% calc ratios 
position = [ 4 5 6 10 11 12 16 17 18 22 23 24];

ratioarray = ones(size(efluxtable,1),length(position))*NaN;

for i= 1:length(position)
    ratioarray(:,i)= eflux_array(:,position(i))./eflux_array(:,(position(i)-3));
end


% pool
ratioarray(ratioarray==0)=NaN; % just to be sure 
ratioarray_pooledMean=ones(size(ratioarray,1),(size(ratioarray,2)/3))*NaN;
ratioarray_pooledSdev=ones(size(ratioarray,1),(size(ratioarray,2)/3))*NaN;

% calc ratio mean 
for col = 0:(size(ratioarray_pooledMean,2)-1)
    ratioarray_pooledMean(:,col+1)=nanmean(ratioarray(:,col*3+1:(col+1)*3),2);
end

% calc ratio sdev
for col = 0:(size(ratioarray_pooledSdev,2)-1)
    ratioarray_pooledSdev(:,col+1)=std(ratioarray(:,col*3+1:(col+1)*3),0,2,'omitnan');
end


% create ratiotable with sdev for tdata
final = ones(size(model.rxnNames,1),(1+size(ratioarray_pooledMean,2)+size(ratioarray_pooledSdev,2)))*NaN;
final = num2cell(final);

final(:,1)= model.rxnNames;
final(:,[2 4 6 8])=num2cell(ratioarray_pooledMean(:,:));
final(:,[3 5 7 9])=num2cell(ratioarray_pooledSdev(:,:));
final=cell2table(final);


% set colnames
varnamez3 = {'atgcode',...
            'RatioT_L_0',...
            'sdev_L_0',...
            'RatioT_L_1',...
            'sdev_L_1',...
            'RatioT_L_7',...
            'sdev_L_7',...
            'RatioT_L_21',...
            'sdev_L_21'};
        
final.Properties.VariableNames=varnamez3;


%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%%%
%%% everything checked by hand, its correct!!! %%%


% write to file 
writetable(final,'final_Tdata_mean_dev.txt','Delimiter','\t');
writetable(final,'final_Tdata_mean_dev.csv','Delimiter','\t');
writetable(final,'final_Tdata_mean_dev.xlsx');


