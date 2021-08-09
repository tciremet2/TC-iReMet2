function [rvCorrStruct] = Store_WTandMT_fluxMatrix (rvCorrStruct, OptResults, cutoff, repetition, kNumbExps)


    %eval (strcat ( 'rvCorrStruct.exponent',num2str(cutoff),'.MT.rep',num2str(repetition),' = OptResults.mutantFluxMatrix;'));
    %eval (strcat ( 'rvCorrStruct.exponent',num2str(cutoff),'.WT.rep',num2str(repetition),' = OptResults.wildtypeFluxMatrix;'));
    
    
    % use main these parts later on as main, since first part only for
    % information 
    for time=0:kNumbExps-1
        
        eval(strcat('rvCorrStruct.e',num2str(cutoff),'.main_WT.t', num2str(time), '(:,', num2str(repetition), ') = OptResults.wildtypeFluxMatrix(:,', num2str(time+1),');'));
        eval(strcat('rvCorrStruct.e',num2str(cutoff),'.main_MT.t', num2str(time), '(:,', num2str(repetition), ') = OptResults.mutantFluxMatrix(:,', num2str(time+1),');'));
        
    end
    
    
end 