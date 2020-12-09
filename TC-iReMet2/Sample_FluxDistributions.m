function [sampledFluxDistributions] = Sample_FluxDistributions( Constraints, optimalValueJ, x0, numberOfRuns )

	
	%%% set important variables
	experiment = 4;
	numberOfRuns = 3; 
		
	
	%%%%%%%%%%%%
	%%% INIT %%%
	%%%%%%%%%%%%
    
	%%% load Results derived from the original optimization
	load Run_1e8.mat % contains all important information from original objective 
	
    %load Run_1e8_all.mat
    %%% set optinput according to experiment as well as number of rxns 
	kNumbRxns = 549;
    eval(strcat('qpOptInput = ResultsPlusConstraints.experiment', num2str(experiment), '.OptInput;'));
	
    % x_L/U
	x_L = qpOptInput.x_L;
	x_U = qpOptInput.x_U;
    
    % name of the new optimization
	name = 'Alpha maximization';
    
    % save old solution
	eval(strcat('init_x_0 = ResultsPlusConstraints.experiment', num2str(experiment), '.Result.x_k;')); %results first init 
	
    % constraint matrix and bounds
    A = qpOptInput.A_total;
	b_L = qpOptInput.b_total_L;
	b_U = qpOptInput.b_total_U;
	
    % save new constraints seperately
	C = A; % C corresponds to the Constraint matrix
	
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ADD ALPHA AND NEW CONSTRAINTS %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	% add col for alpha
    C = [C, zeros(size(C,1),1)]; 
    % add constraint -v_new = -v_old + alpha * d -- -v_new + d*a = -v_old
    newCons = zeros(2*kNumbRxns,size(C,2));
	for i = 1:size(newCons,1)
		newCons(i,i) = -1;
	end
	% add new constraints to overall constraints 
	C = [C; newCons];
	b_L = [b_L; zeros(size(newCons,1),1)]; % for now set to 0
	b_U = [b_U; zeros(size(newCons,1),1)]; % for now set to 0
	
	% add variable for alpha 
	x_L = [x_L; zeros(1)];
    % alternative to fix to find some solution make alpha at least be v sma
    %x_L = [x_L; 0.001];
	x_U = [x_U; 100];
	
    % set bounds for quadratic constraint function - deviation hardcoded change that     
	c_L = 0.99*OptResults.deviations(experiment); % a lot deviation (10%)
	c_U = 1.01*OptResults.deviations(experiment); % a lot deviation (10%)
	
    % set objective function 
	d = zeros(1,size(C,2));
	d(end) = -1; % idx alpha
	
    % init pre while starting vector 
	v_init = [init_x_0; zeros(1)];
        

    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% SAMPLE WITH HIT AND RUN %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	alpha = 0;
	% init pre loop
    v_new = v_init;
	% create flag and counter of iteration
    finished = boolean(0);
    counter = 1;
    % tries to find alpha
    alphaTries=0;
    % create sample matrix
	sampledFluxDistributions = ones(size(v_init,1),numberOfRuns)*NaN;
    % run as often as one specifies number of runs 
    while isequal(finished,0)
       
		% init vold 
		v_old = v_new; % v_new in first iteration is the old optimal solution	
		% set bounds of v_old constraint 
        b_L((size(A,1)+1):end) = -1*v_old(1:2*kNumbRxns); %-0.001;
        b_U((size(A,1)+1):end) = -1*v_old(1:2*kNumbRxns); %+0.001;
                
        tic() 
        while isequal(alpha,0)
            % create random direction vector - direc
            direc = zeros(1,2*kNumbRxns);
            for i = 1:length(direc)
                direc(i) = (rand()-0.5)*1e-3; % create random number between -0.5 and 0.5; does size of direction matter? no i guess 
            end
            % set direction for alpha 
            C((size(A,1)+1):end,end) = direc(:);
            
            % max a
            switch experiment
                case 1
                    Prob = lpconAssign(d, x_L, x_U, name, v_old, C, b_L, b_U, 'Set_SamplingQconstraint_experiment1', [], [], [], c_L, c_U );
                    %Prob = lpconAssign(d, x_L, x_U, name, v_old, C, b_L, b_U, [], [], [], [], [], [] );
                case 2
                    Prob = lpconAssign(d, x_L, x_U, name, v_old, C, b_L, b_U, 'Set_SamplingQconstraint_experiment2', [], [], [], c_L, c_U );
                case 3
                    Prob = lpconAssign(d, x_L, x_U, name, v_old, C, b_L, b_U, 'Set_SamplingQconstraint_experiment3', [], [], [], c_L, c_U );
                case 4
                    Prob = lpconAssign(d, x_L, x_U, name, v_old, C, b_L, b_U, 'Set_SamplingQconstraint_experiment4', [], [], [], c_L, c_U );
            end
            %Prob.SOL.SpecsFile = 'tomlab.specs';
            Result = tomRun('snopt', Prob, 1);
            disp(Result.f_k(end));
            if abs(Result.f_k(end))>1e-4
                alpha = Result.f_k(end);
                disp(alpha);
            end
            % change alpha only if one finds a solution
            %if Result.ExitFlag == 0
            %    alpha = Result.f_k(end);
            %    disp(alpha);
            %end
            alphaTries = alphaTries + 1;
            if alphaTries > 1000
                print('couldn''t find an direction vector pointing into solutionspace after 1000 trials!');
				break
            end
        end
		toc()
        % pick random distance in range of alpha a
        randHowFarMove = 0 + (alpha * rand());
        move = (randHowFarMove * direc)'; 
        move = [move; zeros((size(v_old,1) - (size(move,1)) ),1)];
        % update v
        v_new = v_old + move; 
        % save sample
        sampledFluxDistributions(:,counter) = v_new;
        % set alpha again to 0 for the whileloop condition (depends on how to actually go about it alpha/exitstatus) 
        alpha = 0; 
        
        % finish condition 
        if isequal(counter,numberOfRuns)
            finished = boolean(1);
        end
        
        % add to sample counter 
        counter = counter + 1;
        
	end

end