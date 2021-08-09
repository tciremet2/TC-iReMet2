function [ fluxArrayDayBefore ] = Set_fluxArrayDayBefore (fluxArrayDayBefore, OptResults, day)

    % set WT flux values 
    fluxArrayDayBefore(:,1) = OptResults.wildtypeFluxMatrix(:,day);
    % set MT flux values 
    fluxArrayDayBefore(:,2) = OptResults.mutantFluxMatrix(:,day);

end