function [ OptResults ] = Init_OptResults( kNumbExps, kNumbRxns, kNumbRatioConstr )

%   OptResults.deviations: Euclidian distance between wildtype and mutant
%   OptResults.deviationMatrix: difference flux vectors between wildtype and mutant
%   OptResults.wildtypeFluxMatrix: flux vectors wildtype
%   OptResults.mutantFluxMatrix: flux vectors mutant

    OptResults.deviations              = ones(kNumbExps, 1)*NaN;
    OptResults.deviationMatrix         = ones(kNumbRxns, kNumbExps)*NaN;
    OptResults.wildtypeFluxMatrix      = ones(kNumbRxns, kNumbExps)*NaN;
    OptResults.mutantFluxMatrix        = ones(kNumbRxns, kNumbExps)*NaN;
    OptResults.fluxRangeMatrixWT       = ones(kNumbRxns, 2*kNumbExps)*NaN;
    OptResults.fluxRangeMatrixM        = ones(kNumbRxns, 2*kNumbExps)*NaN;
    OptResults.fluxRangeDiffMatrixWT   = ones(kNumbRxns, kNumbExps)*NaN;
    OptResults.fluxRangeDiffMatrixM    = ones(kNumbRxns, kNumbExps)*NaN;
    
    OptResults.slackValueMatrix        = ones(kNumbRatioConstr, kNumbExps)*NaN;
    OptResults.ratioConstraintMatrix   = ones(kNumbRatioConstr, 2*kNumbExps)*NaN; 

end
