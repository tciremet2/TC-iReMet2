function [ OptInitStruct ] = Init_OptArrays( model, kNumbRxns, kNumbMets )

	Name  = 'Minimization of flux distance';

	% Create matrix defining quadratic form; yielding: (x_i - x_(i+n))^2 + ...
	% + (x_i - x_(i-n))^2 <- first term for i <= n (wildtype entries of combined flux vector), second
	% for i > n (mutant entries)
	F = speye( 2 * kNumbRxns );
	for i = 1 : kNumbRxns 
        F(i, kNumbRxns + i) = -1;
        F(kNumbRxns + i , i) = -1;
	end
	F = 2 * F; % qpconAssign defines the objective  "0.5 * x' * F * x + d' * x"

	d = zeros(2 * kNumbRxns, 1);

	x_0 = zeros(2 * kNumbRxns, 1);

	% x_L = [vmin; vmin];
	% x_U = [vmax; vmax];

	% Linear constraints - steady state flux
	A = [model.S, zeros(kNumbMets, kNumbRxns); zeros(kNumbMets, kNumbRxns), model.S];
	b_L = zeros(2 * kNumbMets, 1);
	b_U = zeros(2 * kNumbMets, 1);
	
    OptInitStruct.Name = Name;
	OptInitStruct.A = A;
	OptInitStruct.b_L = b_L;
	OptInitStruct.b_U = b_U;
    OptInitStruct.d = d;
	OptInitStruct.F = F;
	OptInitStruct.x_0 = x_0;

	% % Linear constraints - biomass constraint
	% A_biomassW = zeros(1, 2*kNumbRxns);
	% A_biomassW(biomass_reaction_id) = 1; % wildtype biomass production
	% b_L_biomassW = fopt;
	% b_U_biomassW = fopt;
	
end