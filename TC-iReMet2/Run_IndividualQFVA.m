function [ fluxRangeWT, fluxRangeM ] = Run_IndividualQFVA( Result, OptResults, ExpCondBoundaries, OptInput, AP, experiment )

    kNumbRxns = size(OptResults.deviationMatrix, 1);
    kNumbRatioConstr = size(OptResults.ratioConstraintMatrix, 1);
    
    examined_rxns = OptResults.qpFVA_examined_rxns;
    
    fluxRangeWT = ones(length(examined_rxns), 2)*NaN;
    fluxRangeM  = ones(length(examined_rxns), 2)*NaN;
    
    x_L = OptInput.x_L;
    x_U = OptInput.x_U;
    
    init_x_0 = Result.x_k;
    x = Result.x_k;
    
    x_L( ExpCondBoundaries.biomRxnIdWT ) = (1-AP.allowed_rel_dev_linear) * Result.x_k( ExpCondBoundaries.biomRxnIdWT );
    x_U( ExpCondBoundaries.biomRxnIdWT ) = (1+AP.allowed_rel_dev_linear) * Result.x_k( ExpCondBoundaries.biomRxnIdWT );
    x_L( kNumbRxns + ExpCondBoundaries.biomRxnIdM ) = (1-AP.allowed_rel_dev_linear) * Result.x_k( kNumbRxns + ExpCondBoundaries.biomRxnIdM );
    x_U( kNumbRxns + ExpCondBoundaries.biomRxnIdM ) = (1+AP.allowed_rel_dev_linear) * Result.x_k( kNumbRxns + ExpCondBoundaries.biomRxnIdM );

    %%%  SLACKS
    if AP.with_slacks == 1
        for slackID = 2*kNumbRxns+1 : 2*kNumbRxns+kNumbRatioConstr 
            if Result.x_k(slackID) >=  0
                x_L(slackID) = (1-AP.allowed_rel_dev_linear)*Result.x_k(slackID);
                x_U(slackID) = (1+AP.allowed_rel_dev_linear)*Result.x_k(slackID);
            else 
                x_L(slackID) = (1+AP.allowed_rel_dev_linear)*Result.x_k(slackID);
                x_U(slackID) = (1-AP.allowed_rel_dev_linear)*Result.x_k(slackID);
            end
        end
    end
    disp('made it after slacks');
    disp(length(x_L));
    
    % Nonlinear constraints
    c_L = (1-AP.allowed_rel_dev_nonlinear) * OptResults.deviations(experiment);
    c_U = (1+AP.allowed_rel_dev_nonlinear) * OptResults.deviations(experiment);

    optDirecs = [1; -1];
    for rxn = 1:(2*length(examined_rxns))
        disp(rxn)
        % Init
        d = zeros(size(OptInput.A_total, 2), 1); % Captures linear parts

        for direc = 1:2
            if rxn <= length(examined_rxns)
                d(examined_rxns(rxn)) = optDirecs(direc);
            else
                d(examined_rxns(rxn-length(examined_rxns))+kNumbRxns) = optDirecs(direc);
            end
            
            %%% CP for everyday case with different functions depending on
            %%% experiment thingy  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            % Create tomlab problem structure - different quadratic constraints since weights of them are time dependend
            switch experiment
                case 1
                    Prob = lpconAssign(d, x_L, x_U, 'Maximization of flux under qp constraint', init_x_0, OptInput.A_total, OptInput.b_total_L, OptInput.b_total_U, 'Set_qpConstraintFunc_experiment1', [], [], [], c_L, c_U );
                case 2
                    Prob = lpconAssign(d, x_L, x_U, 'Maximization of flux under qp constraint', init_x_0, OptInput.A_total, OptInput.b_total_L, OptInput.b_total_U, 'Set_qpConstraintFunc_experiment2', [], [], [], c_L, c_U );
                case 3
                    Prob = lpconAssign(d, x_L, x_U, 'Maximization of flux under qp constraint', init_x_0, OptInput.A_total, OptInput.b_total_L, OptInput.b_total_U, 'Set_qpConstraintFunc_experiment3', [], [], [], c_L, c_U );
                case 4
                    Prob = lpconAssign(d, x_L, x_U, 'Maximization of flux under qp constraint', init_x_0, OptInput.A_total, OptInput.b_total_L, OptInput.b_total_U, 'Set_qpConstraintFunc_experiment4', [], [], [], c_L, c_U );
            end
            
%             % Run a global solver for 10000 function evaluations.
%             Result = tomRun('glcFast', Prob, 1); 

            %%% GLCFAST -- "Since no such constant is used, there is no natural way of defining convergence 
            %%% (except when the optimal function value is known). Therefore glcFast  is run for a predefined 
            %%% number of function evaluations and considers the best
            %%% function value found as the optimal one."
            %%%
            %%% --> the TOMLAB quickguide suggests searching first with
            %%% glcfast (global), then with snopt (local, recognizes
            %%% convergence); however, here, we can already start from a
            %%% known solution -- therefore it makes sense to skip glcfast and use snopt from the beginning!!
%             Prob.x_0 = init_x0; % Use the result from the distance
%             minimization as start vector <- already done!

            % Run SNOPT as a local solver
%             Prob.SOL.SpecsFile = 'tomlab.specs';
            Result = tomRun('snopt', Prob, 1);

            if Result.ExitFlag == 0
                if rxn <= length(examined_rxns)
                    fluxRangeWT(rxn, direc) = Result.f_k;
                else
                    fluxRangeM(rxn-length(examined_rxns), direc) = Result.f_k;
                end
            else
                if rxn <= length(examined_rxns)
                    fluxRangeWT(rxn, direc) = NaN;
                else
                    fluxRangeM(rxn-length(examined_rxns), direc) = NaN;
                end
            end
        end

    end

    fluxRangeWT(:, 2) = -1*fluxRangeWT(:, 2);
    fluxRangeM(:, 2)  = -1*fluxRangeM(:, 2);
    
end
