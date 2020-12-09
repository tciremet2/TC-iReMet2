function [ considered_experiments ] = Determine_ConsideredExperiments( ExpDataSpecific, env_diurnal, experimentname )

    considered_experiments = [];
    
    for experimentID = 1:length(ExpDataSpecific.description)
        % Experiment-specific diurnal conditions
        if env_diurnal == 'day'
            if ~isempty(strfind(ExpDataSpecific.description{2, experimentID}, experimentname)) 
                considered_experiments = [ considered_experiments, experimentID ];
            end
        elseif env_diurnal == 'night'
            if ~isempty(strfind(ExpDataSpecific.description{2, experimentID}, '-20')) % -20 from original approach still there, in own case no 'night data' present
                considered_experiments = [ considered_experiments, experimentID ];
            end
        end
    end
    
end
