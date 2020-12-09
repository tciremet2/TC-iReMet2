function [ Result ] = Run_TomlabWithDifferentStartingValues ( Name, OptInput, kNumbTestX_0 )
 
 % Initialize testing multiple x_0
        found_feasible = 0;
        distance_save = Inf;
        
        % Create tomlab problem structure
        Prob = qpconAssign(OptInput.F, OptInput.d, OptInput.x_L, OptInput.x_U, Name, OptInput.x_0, OptInput.A_total, OptInput.b_total_L, OptInput.b_total_U);
        Prob.SOL.SpecsFile = 'tomlab.specs';
        % Run SNOPT as a local solver
        Result = tomRun('snopt', Prob, 1);
%         Prob = qpAssign(OptInput.F, OptInput.d, OptInput.A_total, OptInput.b_total_L, OptInput.b_total_U, OptInput.x_L, OptInput.x_U);
%         Prob.MIP.cpxControl.QPMETHOD=2;
%         Prob.PriLevOpt = 2;
%         Result = tomRun('cplex',Prob,1);
        
        if Result.ExitFlag == 0
            found_feasible = 1;    
            distance_save = 0.5 * Result.x_k' * OptInput.F * Result.x_k;
            x_0_save = OptInput.x_0;
        end

        if kNumbTestX_0 > 1
            for test_x_0 = 1:kNumbTestX_0
                x_0 = 100 * rand(length(OptInput.x_0), 1);
                % Create tomlab problem structure
                Prob = qpconAssign(OptInput.F, OptInput.d, OptInput.x_L, OptInput.x_U, Name, x_0, OptInput.A_total, OptInput.b_total_L, OptInput.b_total_U);
                Prob.SOL.SpecsFile = 'tomlab.specs';

                % Run SNOPT as a local solver
                Result = tomRun('snopt', Prob, 1);

                if Result.ExitFlag == 0
                    found_feasible = 1;    

                    if 0.5 * Result.x_k' * OptInput.F * Result.x_k < distance_save
                        distance_save = 0.5 * Result.x_k' * OptInput.F * Result.x_k;
                        x_0_save = x_0;
                    end
                end
            end

            if found_feasible == 1
                x_0 = x_0_save;
                Prob = qpconAssign(OptInput.F, OptInput.d, OptInput.x_L, OptInput.x_U, Name, x_0, OptInput.A_total, OptInput.b_total_L, OptInput.b_total_U);
                Prob.SOL.SpecsFile = 'tomlab.specs';
                % Run SNOPT as a local solver
                Result = tomRun('snopt', Prob, 1);
            end 
        end
        
end