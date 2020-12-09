function [model] = Correct_model_mistakes(model)

    over9charlengthIndices = ~cellfun(@isempty,regexp(model.genes,'\w\w\w\w\w\w\w\w\w\w'));
    mistakes_modelgenes = model.genes(over9charlengthIndices); % --> 2 occurences: one rule, one concatenated names(count2)
    % rule: single atgs not present therefore add them under the assumption it is just a writing mistake; AT-codes part of atp-citrate lyase
    model.genes(strcmp(model.genes,mistakes_modelgenes{2}))={'-'}; % exclude name  %1
    % AT3G06650 is not present --> add it to names
    model.genes{end+1}='AT3G06650'; %1 
    % AT5G49460 is not present --> add it to names
    model.genes{end+1}='AT5G49460'; %1
    % 2 concatenated names: single atgs present therefore just eliminate the rule 
    model.genes(strcmp(model.genes(:,1),'AT1G31230AT4G19710'))={'-'}; % allowed because there is alreadz a '-' in the vector %2
    %%% model.grRules - problems found when doing efluxlike method und hardcoded back  
    model.grRules(130)= {'4*(AT1G10670 OR AT1G60810 OR AT1G09430) AND 4*(AT3G06650 OR AT5G49460)'}; % writing mistake i guess %1
    % no operator in rule, in araCore publication published like this as well
    % therefore checked databases and it seems to be hexamer
    % thats why OR operator gets substituted;
    % ask for: reason behind it?
    model.grRules(235)= {'4*(AT1G31230 OR AT4G19710)'}; % in brenda homohexamer --> OR        %2 atgs are the same as in the concatenated name
    model.grRules(236)= {'4*(AT1G31230 OR AT4G19710)'}; % fuer 2te atg nichts vorhanden       %2 atgs are the same as in the concatenated name

end