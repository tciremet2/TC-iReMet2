function [ distance ] = Set_qpConstraintFunc( x, Prob )

    kNumbRxns = 549; % no possibility to assign this in a different way

    % Create matrix defining quadratic form; yielding: (x_i - x_(i+n))^2 + ...
    % + (x_i - x_(i-n))^2 <- first term for i <= n (wildtype entries of combined flux vector), second
    % for i > n (mutant entries)
    F = speye( 2 * kNumbRxns );
    for i = 1 : kNumbRxns 
      F(i, kNumbRxns + i) = -1;
      F(kNumbRxns + i , i) = -1;
    end

    distance = x(1:(2*kNumbRxns))' * F * x(1:(2*kNumbRxns)); 
    

end 