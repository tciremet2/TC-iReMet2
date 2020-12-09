function [ distance ] = Set_qpConstraintFunc_experiment2( x, Prob )

    % set experimentnumber
    exp = 2;
    
    %%% initialize and read in necessary information 
    weighttable = Read_weightTable();
    weights = weighttable(:,exp);
    % CHANGE THIS ACCORDING TO USED CUTOFF!!!!!!!!!!!!!!!!!!!
    %linear_constraintPart = dlmread('/home/mpimp-golm.mpg.de/pries1347/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/day/without_constCofactorRatios/with_slacks/with_slack_minimization/1e8/rep1/linear_objectivePart_experiments.dat');
    %const_add = dlmread('/home/mpimp-golm.mpg.de/pries1347/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/day/without_constCofactorRatios/with_slacks/with_slack_minimization/1e8/rep1/const_add_objective.dat');
    kNumbRxns = 549; % no possibility to assign this in a different way
    
    load('optintputest.mat')

    % Create matrix defining quadratic form; yielding: (x_i - x_(i+n))^2 + ...
    % + (x_i - x_(i-n))^2 <- first term for i <= n (wildtype entries of combined flux vector), second
    % for i > n (mutant entries)
    F = speye( 2 * kNumbRxns );
    for i = 1 : kNumbRxns 
      F(i, kNumbRxns + i) = -1;
      F(kNumbRxns + i , i) = -1;
    end
    
    % calculate weights for day one
    %%% Init weights and constants 
    % set weights: alpha - a, beta - b, gamma - g 
    a = weights(1);
    b = weights(2);
    g = weights(3);
    % alpha+beta and alpha+gamma for quadratic
    ab = a+b;
    ag = a+g;

    % adjust the objective function according to weights of that experiemnt
    for i=1:kNumbRxns
           F(i,i) = F(i,i)*ab;
           F(kNumbRxns+i,kNumbRxns+i) = F(kNumbRxns+i,kNumbRxns+i)*ag;
           F(i,kNumbRxns+i) = F(i,kNumbRxns+i)*a;
           F(kNumbRxns+i,i) = F(kNumbRxns+i,i)*a;
    end
    
    % add slackparts to F to match dimensions
    F = [F, zeros(size(F,1),length(x)-2*kNumbRxns)];
    F = [F; zeros(length(x)-2*kNumbRxns,size(F,2))];
        
    % set linear part of the quadratic constraint derived from the original
    % objective function
    %d = linear_constraintPart(:,exp);
    
    % calculate distance 
    distance = x' * F * x + OptInput.d' * x + OptInput.const_add;
    
end