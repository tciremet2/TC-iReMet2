% CP concat table for leaf transcript data 
% here geneannotation and datapoints get put together for further analysis
% via eflux and ratio calculations 

%%% read in tables

% import leaf3 via right click
%
% IMPORTANT: MARK 3RD COLUMN AS TEXT NOT NUMERIC!!!!!!
%
%%%%%  IMPORT RIGHT CLICK %%%%
 
%%% use transcriptleaf3 as data
data = transcriptleaftable3;
dat  = table2array(data);

%%% remove transcript 
% remove .x of atgcodes 
tmp(:) = strtok(dat(:,3),'.');
tmp=tmp';
dat(:,3)=tmp;

%%% create mapping and prestuff for aggregation 
[atgs, ~, subs] = unique(dat(:,3), 'rows');
% create map 
mapping = cell(length(dat(:,3)),2);
% fill first col with atg codes
for i=1:length(dat(:,3))
    mapping{i,1}=dat(i,3);
end
% fill second col with subs
for i=1:length(subs)
    mapping{i,2}=subs(i);
end


%%% convert data to numeric array so it can be summarized
dat_ID = [subs,dat];
dat_ID = array2table(dat_ID);
numericdat = table2array(dat_ID(:,[1 4 7:end]));
test = table2array(dat_ID(:,[1 4 7:end]));
test(:,8) = '0'; % replace missing with 0 since they dont normally occur
test(:,2) = []; % remove col with strings that cant be converted to double
A = array2table(test);
A = table2array(A);
A = str2double(A);
A = num2cell(A);
B = cell2mat(A);


% keep on aggregation, vector length of data col  
smrz = 1:size(B,2);

% aggregation 
aggregatetable = zeros(size(atgs, 1), numel(smrz));
for colidx = 1:numel(smrz)
   aggregatetable(:, colidx) = accumarray(subs, B(:, smrz(colidx)), [], @mean);
end

final = aggregatetable;
% final is corrected, checked with examples picked by hand 

%%% convert final back to table that can be saved 
test = num2cell(final);
test2 = test;

% only run this if you have time :^), exchanges numbers with atgcodes 
for aggreID=1:size(test,1)
    for mapID=1:size(mapping,1)
        if (isequal(test(aggreID,1),mapping(mapID,2)))
            test2(aggreID,1)=mapping(mapID,1);
        end
    end
end

% set header 
varnamez2 = {'atgcode',...
'Col0_L_0_1','Col0_L_0_2','Col0_L_0_3',...
'reil11reil21_L_0_1','reil11reil21_L_0_2','reil11reil21_L_0_3',...
'Col0_L_1_1','Col0_L_1_2','Col0_L_1_3',...
'reil11reil21_L_1_1','reil11reil21_L_1_2','reil11reil21_L_1_3',...
'Col0_L_7_1','Col0_L_7_2','Col0_L_7_3',...
'reil11reil21_L_7_1','reil11reil21_L_7_2','reil11reil21_L_7_3',...
'Col0_L_21_1','Col0_L_21_2','Col0_L_21_3',...
'reil11reil21_L_21_1','reil11reil21_L_21_2','reil11reil21_L_21_3'};

% create final table 
ultfinal = test2;
ultfinal2 = cell2table(ultfinal,'VariableNames',varnamez2);

% write final table, the 0 col has to exchanged by NAs 
writetable(ultfinal2,'final.csv');
writetable(ultfinal2,'final.txt');
writetable(ultfinal2,'final.xlsx');

% average replicates
nextA = table2array(ultfinal2(:,2:end));
nextA(nextA==0)=NaN;
A_pooledMean=zeros(size(ultfinal2,1),((size(ultfinal2,2)-1)/3));
for col = 0:((size(nextA,2)/3)-1)
    A_pooledMean(:,col+1)=nanmean(nextA(:,col*3+1:(col+1)*3),2);
end

A_pooledMean=array2table(A_pooledMean);

oi = [ultfinal2(:,1),A_pooledMean];

varnamez3 = {'atgcode',...
            'Col0_L_0',...
            'reil11reil21_L_0',...
            'Col0_L_1',...
            'reil11reil21_L_1',...
            'Col0_L_7',...
            'reil11reil21_L_7',...
            'Col0_L_21',...
            'reil11reil21_L_21'};
    
%sd
oi.Properties.VariableNames=varnamez3;


writetable(oi,'final_pooled.csv');
writetable(oi,'final_pooled.txt');
writetable(oi,'final_pooled.xlsx');
    

%{
%%% create ratios here 
% WT/MT or MT/WT --> you have to set this by coding it - not nice but works
% init table to save ratios in 
ratioTable = zeros(size(oi,1),(size(oi,2)-1)/2);
calcmat = table2array(oi(:,2:end));
% calc ratios
for col = (0:(size(calcmat,2)/2)-1)
    ratioTable(:,col+1)=calcmat(:,col*2+1)./calcmat(:,(col+1)*2);
end
% make ratio table 
tasty = [oi(:,1),array2table(ratioTable)];
% set varnamez
tasty.Properties.VariableNames={'ATGcode','ratio_L_0','ratio_L_1','ratio_L_7','ratio_L_21'};


writetable(tasty,'transcript_ratios.csv');
writetable(tasty,'transcript_ratios.txt');
writetable(tasty,'transcript_ratios.xlsx');

%}



    
    