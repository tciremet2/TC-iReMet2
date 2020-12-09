function [rvCorrStruct] = Init_rvCorrStruct(maxratio_array, AnalysisParameters, kNumbRxns, kNumbExps )

    % create a matrix for correlation of fluxes for every exponents reruns
    for exponentOfCutoff = 1:length(maxratio_array)
        
        for reps = 1:AnalysisParameters.kNumbRepetition

           eval(strcat('rvCorrStruct.exponent',num2str(exponentOfCutoff),'.WT.rep',num2str(reps),' = ones(kNumbRxns,kNumbExps)*NaN;'));
           eval(strcat('rvCorrStruct.exponent',num2str(exponentOfCutoff),'.MT.rep',num2str(reps),' = ones(kNumbRxns,kNumbExps)*NaN;'));
           
        end
    end
    
    % init hauptteil 
    for exponentcutoff = 1:length(maxratio_array)
        for day=0:kNumbExps-1

           eval(strcat('rvCorrStruct.e',num2str(exponentcutoff),'.main_WT.t',num2str(day),'=ones(kNumbRxns,AnalysisParameters.kNumbRepetition)*NaN;' ));
           eval(strcat('rvCorrStruct.e',num2str(exponentcutoff),'.main_MT.t',num2str(day),'=ones(kNumbRxns,AnalysisParameters.kNumbRepetition)*NaN;' ));

        end
    end
    

end 