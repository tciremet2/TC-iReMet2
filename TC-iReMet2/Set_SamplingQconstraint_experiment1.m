function [distance] = Set_SamplingQconstraint_experiment1(x, Prob)


	% set experimentnumber
    exp = 1;
	
	%%% initialize and read in necessary information 
    weighttable = Read_weightTable();
    weights = weighttable(:,exp);
    kNumbRxns = 549;
	
	% load inputvariables
	%load Set_SamplingQconstraint_experiment1.mat
    load Set_SamplingQconstraint_experiment1_2.mat
	
	% Create matrix defining quadratic form; yielding: (x_i - x_(i+n))^2 + ...
    % + (x_i - x_(i-n))^2 <- first term for i <= n (wildtype entries of combined flux vector), second
    % for i > n (mutant entries)
    F = speye( 2 * kNumbRxns );
    for i = 1 : kNumbRxns 
      F(i, kNumbRxns + i) = -1;
      F(kNumbRxns + i , i) = -1;
    end
	
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
    
    % add slackparts to F to match dimensions (also accounts for alpha)
    F = [F, zeros(size(F,1),length(x)-2*kNumbRxns)];
    F = [F; zeros(length(x)-2*kNumbRxns,size(F,2))];
	
	% calculate distance 
    distance = x' * F * x + [OptInput1.d;0]' * x + OptInput1.const_add;
	
	
end