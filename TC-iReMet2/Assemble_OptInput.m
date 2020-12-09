function [ OptInput ] = Assemble_OptInput ( AnalysisParameters, OptInit, kNumbRatioConstr, ...
		ExpCondBoundaries, RatioConstraint, kNumbRxns, fluxArrayDayBefore, day, weightsObjectiveDayspecific )
    
    F   = OptInit.F;
    A   = OptInit.A;
    d   = OptInit.d;
    x_0 = OptInit.x_0;
    
    if AnalysisParameters.with_slacks == 1

        % Adding slack variables to test elasticity needed to turn
        % infeasible linear constraints into feasible
        %
        % 2 slack variables for each ratio constraint: one for the ratio
        % constraint and one additional ones to bound the slack variable, in
        % order to minimize sum of absolutes of ratio constraint slacks

        % Extend F
        F = [ F, zeros(size(F, 1), AnalysisParameters.kNumbSlacksPerRatioConstr * kNumbRatioConstr) ]; 
        F = [ F; zeros(AnalysisParameters.kNumbSlacksPerRatioConstr * kNumbRatioConstr, size(F, 2)) ];

        % Extend A  
        A = [ A, zeros(size(A, 1), AnalysisParameters.kNumbSlacksPerRatioConstr * kNumbRatioConstr) ];

        % Extend d    
        if AnalysisParameters.with_slack_minimization == 1
            %d = [d; zeros(kNumbRatioConstr, 1); AnalysisParameters.slack_weighting * ones( kNumbRatioConstr, 1 )]; 
            % CP upper version original ( change for slack s`)
            d = [d; AnalysisParameters.slack_weighting * ones( kNumbRatioConstr, 1 )]; 
        else
            d = [d; zeros(kNumbRatioConstr, 1)]; 
        end

        % Extend x_0
        x_0 = [x_0; zeros( AnalysisParameters.kNumbSlacksPerRatioConstr * kNumbRatioConstr, 1 )];

        % Slack boundaries
        if AnalysisParameters.with_slack_minimization == 1
            %vmin_slack = [ones(kNumbRatioConstr, 1) * -Inf; zeros( kNumbRatioConstr, 1 )];
            % CP upper part original version 
            vmin_slack = zeros( kNumbRatioConstr, 1 );
        else
            %vmin_slack = [ones(kNumbRatioConstr, 1) * -Inf];
            % CP upper part original version (only needs one slack variable) 
            vmin_slack = zeros(kNumbRatioConstr, 1);
        end
        vmax_slack = ones(AnalysisParameters.kNumbSlacksPerRatioConstr * kNumbRatioConstr, 1) * Inf;

        A_total = [A; RatioConstraint.ratio_constraint_matrix];
        b_total_L = [OptInit.b_L; RatioConstraint.ratio_b_L];
        b_total_U = [OptInit.b_U; RatioConstraint.ratio_b_U];

        x_L = [ExpCondBoundaries.vminWT; ExpCondBoundaries.vminM; vmin_slack];
        x_U = [ExpCondBoundaries.vmaxWT; ExpCondBoundaries.vmaxM; vmax_slack];

    else

        A_total = [A; RatioConstraint.ratio_constraint_matrix];
        b_total_L = [OptInit.b_L; RatioConstraint.ratio_b_L];
        b_total_U = [OptInit.b_U; RatioConstraint.ratio_b_U];

        x_L = [ExpCondBoundaries.vminWT; ExpCondBoundaries.vminM];
        x_U = [ExpCondBoundaries.vmaxWT; ExpCondBoundaries.vmaxM];

    end

        
    if isequal(AnalysisParameters.with_daySpecificObjective,1)
        
       %%% Init weights and constants 
       % set weights: alpha - a, beta - b, gamma - g 
       a = weightsObjectiveDayspecific(1,day);
       b = weightsObjectiveDayspecific(2,day);
       g = weightsObjectiveDayspecific(3,day);
       % alpha+beta and alpha+gamma for quadratic
       ab = a+b;
       ag = a+g;              

       %%% Set F according to new weights  
       for i=1:kNumbRxns
           F(i,i) = F(i,i)*ab;
           F(kNumbRxns+i,kNumbRxns+i) = F(kNumbRxns+i,kNumbRxns+i)*ag;
           F(i,kNumbRxns+i) = F(i,kNumbRxns+i)*a;
           F(kNumbRxns+i,i) = F(kNumbRxns+i,i)*a;
       end
       
       %%% set d
       kWT = fluxArrayDayBefore(:,1);
       kMT = fluxArrayDayBefore(:,2);
       % derp 
       kWTnew = (-b*2)*kWT;
       kMTnew = (-g*2)*kMT;
       % set d
       d = [kWTnew;kMTnew;AnalysisParameters.slack_weighting * ones( kNumbRatioConstr, 1 )];
       
       % calculate constant addition 
       const_add = sum(b*(kWT.^2)) + sum(g*(kMT.^2));
       
    end   
    
    
    
    OptInput.A_total = A_total;
    OptInput.b_total_L = b_total_L;
    OptInput.b_total_U = b_total_U;
    OptInput.x_L = x_L;
    OptInput.x_U = x_U;
    OptInput.F = F;
    OptInput.A = A;
    OptInput.d = d;
    OptInput.x_0 = x_0;
    OptInput.const_add = const_add;
    
    
    
        
end