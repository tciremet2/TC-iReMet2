function [ OptResults ] = Store_OptResults( OptResults, Result, OptInput, experimentID, kNumbRxns, AnalysisParameters, ConstraintInfo, RatioConstraint )

        if Result.ExitFlag == 0
            OptResults.deviations(experimentID)             = 0.5 * Result.x_k' * OptInput.F * Result.x_k + OptInput.d'  * Result.x_k + OptInput.const_add;
            OptResults.wildtypeFluxMatrix(:, experimentID)  = Result.x_k(1:kNumbRxns);
            OptResults.mutantFluxMatrix(:, experimentID)    = Result.x_k((kNumbRxns+1):(2*kNumbRxns));
            OptResults.deviationMatrix(:, experimentID)     = OptResults.wildtypeFluxMatrix(:, experimentID) - OptResults.mutantFluxMatrix(:, experimentID);
            
            if AnalysisParameters.with_slacks == 1
                OptResults.slackValueMatrix(:, experimentID) = Result.x_k(2*kNumbRxns+1 : 2*kNumbRxns+ConstraintInfo.kNumbRatioConstr);
            end

            OptResults.ratioConstraintMatrix(:, [(2*experimentID-1), 2*experimentID]) = RatioConstraint.constraint_matrix_for_readout;
            
        end

end
