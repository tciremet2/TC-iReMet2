% CP script to calculate ratios mean and sdev

clear;
clc;


%%% FIRST TABLE 
% read in table
ratio_e1 = readtable('MTWT_ratio_mets_0d1d7d21d_e1.csv');
% save names
metnames = ratio_e1(:,1);
summarynames = {'mt_wt_0d_e1_meanratio','mt_wt_0d_e1_sdev',...
                'mt_wt_1d_e1_meanratio','mt_wt_1d_e1_sdev',...
                'mt_wt_7d_e1_meanratio','mt_wt_7d_e1_sdev',...
                'mt_wt_21d_e1_meanratio','mt_wt_21d_e1_sdev'};

% convert and remove anomalies
ratio_e1 = table2cell(ratio_e1);
ratio_e1(strcmp(ratio_e1,'#DIV/0!'))= {NaN};
ratio_e1(strcmp(ratio_e1,'0'))= {NaN};

% make number array to be able to make calculation 
ratio_e1 = str2double(ratio_e1);
ratio_e1 = ratio_e1(:,2:end);

% calculate mean and sdev for all ratios 
x=zeros(size(ratio_e1,1),size(ratio_e1,2))*NaN;
for row = 1:size(ratio_e1,1)
    for colclust = 1:3:size(ratio_e1,2)
        x(row,colclust)=nanmean(ratio_e1(row,(colclust:colclust+2)));
        x(row,(colclust+1))=std(ratio_e1(row,(colclust:colclust+2)),'omitnan');
    end
end
% remove unnecessary column
x=x(:,[1 2 4 5 7 8 10 11]);
% name table correctly
ratio_e1_mean_sdev=[metnames,array2table(x)];
ratio_e1_mean_sdev.Properties.VariableNames(2:end)=summarynames;

%1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from the dictionary table as well 
dic = readtable('model_data_dictionary_corrected_processedByAlex_bigVersion.xls');
ratio_e1_mean_sdev = ratio_e1_mean_sdev(2:end,:);

% are measured things the same
testi = table2cell(ratio_e1_mean_sdev(:,1))
testi2 = table2cell(dic(:,1))
i = sort(testi);
j = sort(testi2);
% sort cell table according to stufferinos 
testi = table2cell(ratio_e1_mean_sdev)
testi2 = table2cell(dic)
[~, ix] = sort(testi(:,1));   % sort the first column from 2nd row onwards and get the indices
testi(:,:) = testi(ix,:)
[~, ix] = sort(testi2(:,1))   % sort the first column from 2nd row onwards and get the indices
testi2(:,:) = testi2(ix,:)
% find idx of removers 
test=cell2mat(testi(:,2:end));
x = (1:size(test,1))';
test=[x,test];
test(any(isnan(test), 2), :) = [];
% save idx of left over rows 
leftoveridx = test(:,1);
% ratio table 
testratio = testi;
testratio = testi(leftoveridx,:);
testratio = cell2table(testratio);
% dic table 
testdic = testi2;
testdic = testi2(leftoveridx,:);
testdic = cell2table(testdic);
% write tables 
writetable(testratio,'MTWT_ratio_mets_e1_mean_sdev_sorted.txt');
writetable(testdic,'MTWT_ratio_mets_dictionary_e1_sorted.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




















%%% 2nd TABLE 
% read in table
ratio_e2 = readtable('MTWT_ratio_mets_0d1d7d21d_e2.csv');
% save names
metnames = ratio_e2(:,1);
summarynames = {'mt_wt_0d_e1_meanratio','mt_wt_0d_e1_sdev',...
                'mt_wt_1d_e1_meanratio','mt_wt_1d_e1_sdev',...
                'mt_wt_7d_e1_meanratio','mt_wt_7d_e1_sdev',...
                'mt_wt_21d_e1_meanratio','mt_wt_21d_e1_sdev'};

% convert and remove anomalies
ratio_e2 = table2cell(ratio_e2);
ratio_e2(strcmp(ratio_e2,'#DIV/0!'))= {NaN};
ratio_e2(strcmp(ratio_e2,'0'))= {NaN};

% make number array to be able to make calculation 
ratio_e2 = str2double(ratio_e2);
ratio_e2 = ratio_e2(:,2:end);

% calculate mean and sdev for all ratios 
x=zeros(size(ratio_e2,1),size(ratio_e2,2))*NaN;
for row = 1:size(ratio_e2,1)
    for colclust = 1:3:size(ratio_e2,2)
        x(row,colclust)=nanmean(ratio_e2(row,(colclust:colclust+2)));
        x(row,(colclust+1))=std(ratio_e2(row,(colclust:colclust+2)),'omitnan');
    end
end
% remove ueberfluessige column
x=x(:,[1 2 4 5 7 8 10 11]);
% name table correctly
ratio_e2_mean_sdev=[metnames,array2table(x)];
ratio_e2_mean_sdev.Properties.VariableNames(2:end)=summarynames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from the dictionary table as well 
dic = readtable('model_data_dictionary_corrected_processedByAlex_bigVersion.xls');
ratio_e2_mean_sdev = ratio_e2_mean_sdev(2:end,:);

% are measured things the same
testi = table2cell(ratio_e2_mean_sdev(:,1))
testi2 = table2cell(dic(:,1))
i = sort(testi);
j = sort(testi2);
% sort cell table according to stufferinos 
testi = table2cell(ratio_e2_mean_sdev)
testi2 = table2cell(dic)
[~, ix] = sort(testi(:,1));   % sort the first column from 2nd row onwards and get the indices
testi(:,:) = testi(ix,:)
[~, ix] = sort(testi2(:,1))   % sort the first column from 2nd row onwards and get the indices
testi2(:,:) = testi2(ix,:)
% find idx of removers 
test=cell2mat(testi(:,2:end));
x = (1:size(test,1))';
test=[x,test];
test(any(isnan(test), 2), :) = [];
% save idx of left over rows 
leftoveridx = test(:,1);
% ratio table 
testratio = testi;
testratio = testi(leftoveridx,:);
testratio = cell2table(testratio);
% dic table 
testdic = testi2;
testdic = testi2(leftoveridx,:);
testdic = cell2table(testdic);
% write tables 
writetable(testratio,'MTWT_ratio_mets_e2_mean_sdev_sorted.txt');
writetable(testdic,'MTWT_ratio_mets_dictionary_e2_sorted.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%























%%% THIRD TABLE 
% read in table
ratio_e1e2 = readtable('MTWT_ratio_mets_0d1d7d21d_e1e2.csv');
% save names
metnames = ratio_e1e2(:,1);
summarynames = {'mt_wt_0d_e1e2_meanratio','mt_wt_0d_e1e2_sdev',...
                'mt_wt_1d_e1e2_meanratio','mt_wt_1d_e1e2_sdev',...
                'mt_wt_7d_e1e2_meanratio','mt_wt_7d_e1e2_sdev',...
                'mt_wt_21d_e1e2_meanratio','mt_wt_21d_e1e2_sdev'};

% convert and remove anomalies
ratio_e1e2 = table2cell(ratio_e1e2);
ratio_e1e2(strcmp(ratio_e1e2,'#DIV/0!'))= {NaN};
ratio_e1e2(strcmp(ratio_e1e2,'0'))= {NaN};
% make number array to be able to make calculation 
ratio_e1e2 = str2double(ratio_e1e2);
ratio_e1e2 = ratio_e1e2(:,2:end);

% calculate mean and sdev for all ratios 
x=zeros(size(ratio_e1e2,1),size(ratio_e1e2,2))*NaN;
for row = 1:size(ratio_e1e2,1)
    for colclust = 1:6:size(ratio_e1e2,2)
        x(row,colclust)=nanmean(ratio_e1e2(row,(colclust:colclust+5)));
        x(row,(colclust+1))=std(ratio_e1e2(row,(colclust:colclust+5)),'omitnan');
    end
end
% remove ueberfluessige column
x=x(:,[1 2 7 8 13 14 19 20]);
% name table correctly
ratio_e1e2_mean_sdev=[metnames,array2table(x)];
ratio_e1e2_mean_sdev.Properties.VariableNames(2:end)=summarynames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from the dictionary table as well 
dic = readtable('model_data_dictionary_corrected_processedByAlex_bigVersion.xls');
ratio_e1e2_mean_sdev = ratio_e1e2_mean_sdev(2:end,:);

% are measured things the same
testi = table2cell(ratio_e1e2_mean_sdev(:,1))
testi2 = table2cell(dic(:,1))
i = sort(testi);
j = sort(testi2);
% sort cell table according to stufferinos 
testi = table2cell(ratio_e1e2_mean_sdev)
testi2 = table2cell(dic)
[~, ix] = sort(testi(:,1));   % sort the first column from 2nd row onwards and get the indices
testi(:,:) = testi(ix,:)
[~, ix] = sort(testi2(:,1))   % sort the first column from 2nd row onwards and get the indices
testi2(:,:) = testi2(ix,:)
% find idx of removers 
test=cell2mat(testi(:,2:end));
x = (1:size(test,1))';
test=[x,test];
test(any(isnan(test), 2), :) = [];
% save idx of left over rows 
leftoveridx = test(:,1);
% ratio table 
testratio = testi;
testratio = testi(leftoveridx,:);
testratio = cell2table(testratio);
% dic table 
testdic = testi2;
testdic = testi2(leftoveridx,:);
testdic = cell2table(testdic);
% write tables 
writetable(testratio,'MTWT_ratio_mets_e1e2_mean_sdev_sorted.txt');
writetable(testdic,'MTWT_ratio_mets_dictionary_e1e2_sorted.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
