function [] = Write_Additionalinformation(Information, AnalysisParameters,resultPath,cutoff)


% header names
pattern = {'d0_kNumbMadeConstraints','d0_usedminratio','d0_usedmaxratio','d0_allratios',...
            'd1_kNumbMadeConstraints','d1_usedminratio','d1_usedmaxratio','d1_allratios',...
            'd7_kNumbMadeConstraints','d7_usedminratio','d7_usedmaxratio','d7_allratios',...
            'd21_kNumbMadeConstraints','d21_usedminratio','d21_usedmaxratio','d21_allratios'};
        
% init final table, 1000 arbitrarily chosen, change(hardcoded)
final = ones(10000,4*4)*NaN;

%%%
% counts loop interation
counter = 1;
% blocks of 4 
for i=1:4:16
    
    final(1,i) = Information.numberofconstraintsmade(1,counter);
    final(1:length(Information.usedminratio),i+1) = Information.usedminratio(:,counter);
    final(1:length(Information.usedmaxratio),i+2) = Information.usedmaxratio(:,counter);
    final(1:length(Information.totalmaxratio),i+3) = Information.totalmaxratio(:,counter);
    
    counter = counter + 1;
end

% remove empty rows
final(all(isnan(final),2),:) = [];
% create header and data
new = num2cell(final);
new = [pattern; new];
final = new;
% convert data to a writeable format and write to file 
T = cell2table(final(2:end,:),'VariableNames',final(1,:));
%writetable(T,strcat(resultPath, 'Additional_Information_', AnalysisParameters.experimentname,'_',AnalysisParameters.maxratiostr,...
%                    '_',num2str(AnalysisParameters.min_mutant_BM_prod), '.csv'), 'Delimiter','\t');

if isequaln(AnalysisParameters.with_rangeBiom, 1)
    
    % fraction table read in
    dummy_biompercent_time = table2array(readtable('/home/mpimp-golm.mpg.de/pries1347/winhome/main/MAIN-uebergabe/TC-iReMet2/Data/percentageOfBiomassProductionAtEachTimepoint2.csv'));
      
    
    if isequaln(AnalysisParameters.devrangeBiom, 0 )
        
        writetable(T,strcat(resultPath, 'Additional_Information_', AnalysisParameters.experimentname,'_','1e',num2str(cutoff),...
                    '_',strcat(num2str(dummy_biompercent_time(1,1)),'bm',num2str(dummy_biompercent_time(1,2)),'bm',...
                    num2str(dummy_biompercent_time(1,3)),'bm',num2str(dummy_biompercent_time(1,4))), '.csv'), 'Delimiter','\t');
        
        
    else 
                    
                
        writetable(T,strcat(resultPath, 'Additional_Information_', AnalysisParameters.experimentname,'_','1e',num2str(cutoff),...
                    '_',strcat(num2str(dummy_biompercent_time(1,1)),'bm',num2str(dummy_biompercent_time(1,2)),'bm',...
                    num2str(dummy_biompercent_time(1,3)),'bm',num2str(dummy_biompercent_time(1,4)),...
                    '_range_',num2str(AnalysisParameters.devrangeBiom)), '.csv'), 'Delimiter','\t');         
                
                
    end

else
    
    writetable(T,strcat(resultPath, 'Additional_Information_', AnalysisParameters.experimentname,'_','1e',num2str(cutoff),...
                    '_',num2str(AnalysisParameters.min_mutant_BM_prod), '.csv'), 'Delimiter','\t');
    
    
end

end