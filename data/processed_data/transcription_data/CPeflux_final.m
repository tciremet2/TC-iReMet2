function vmax_eflux = CPeflux_final(model,modellogic,adjustdata,normalize,keepmodel,tdata,excludenonmeasured)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Christopher Pries                                                   %%%    
%%% 788054 University of Potsdam                                        %%%
%%% 30.11.2018                                                          %%%
%%%---------------------------------------------------------------------%%%
%%% iReMetFlux                                                          %%%
%%% CPeflux_final - calculations of enzymeconcentration based on eflux  %%%
%%% For further information about the approach and why its used:        %%%
%%% journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000489% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %%%
% INPUT:                                                                %%%
% model - GEM model, with genenames and grRules                         %%%
% modellogic - modellogic=1 => A OR B, A+B                              %%%
%              modellogic=0 => A OR B, max(A,B)                         %%%
% adjustdata - logical 1 or 0, true if rules are noted as AND/OR        %%%
% normalize  - logical 1 or 0, true results in normalization of the     %%%
%              final values to the max value derived from the eflux     %%%
%              calculations                                             %%%
% keepmodel  - logical 1 or 0, true results in substitution of non      %%%
%              solved rules ub to their original ub                     %%%
% tdata - transcription data in the format of 'atgcode,                 %%%
%         measurement(to do eflux calculations on)'                     %%%
% excludenonmeasured - logical 1 or 0, true if one wants to exclude     %%%
%                      unmeasured at-codes from the eflux calculations  %%%
%-----------------------------------------------------------------------%%%
% OUTPUT:                                                               %%%
% vmax_eflux - vector containing the upper bounds calculated by the     %%% 
%              eflux approach based on the model and dataargument       %%%
%-----------------------------------------------------------------------%%%
% PROCEDURE:                                                            %%%
% 1) Sets variables from model                                          %%%
% 2) Adjust the rules based on the notation                             %%%
% 3) Create a mapping of all genes in the model to their corresponding  %%%
%    measured values, if a gene was not measured, the original name     %%%
%    from the gene names was kept ( important because in 4th step       %%% 
%    at-codes are substituted by their corresponding values based on    %%%
%    the mapping from 3).                                               %%%
% 4) Exchange genenames in the rules with their corresponding values    %%%
%    based on the mapping in 3).                                        %%%
% 5) Store rules on seperate object, to save 'rules'                    %%%
% 6) Remove non measured atcodes from the rules and calculations        %%%
% 7) Resolve and calculate rules                                        %%%
% 8) Normalize calculated rules to the maximum value                    %%%
% 9) Set unresolved rules to their original ub value if specified       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%
%%% 1 %%%
%%%%%%%%%
%%% 1) initializes important variables 
colspec = 2; % 2 because the numbers for calculation are in the 2nd col 
vmax = model.ub;
rules = model.grRules;
%{ 
comments of variables used to make handling easier
%modellogic = 1;
%adjustdata = 1;
%normalize = 0;
%tdata = test;
%excludenonmeasured = 1;
%}


%%%%%%%%%
%%% 2 %%%
%%%%%%%%%
%%% 2) Exchange the notation of operators based on adjustdata  
if adjustdata==1
    for ruleID=1:length(rules)
        rules(ruleID)=strrep(rules(ruleID),'AND','&');
        rules(ruleID)=strrep(rules(ruleID),'OR','|');
    end
end

%%%%%%%%%
%%% 3 %%%
%%%%%%%%%
%%% 3) Create a mapping of modelgenes to the corresponding values in tdata 
% set modelgenenames, measuredgenenames and their corresponding values 
modelgenes = model.genes;
datanames  = tdata(:,1);
datavalues = tdata(:,colspec);
% create mapping, where the 1st columns corresponds to the modelnames, the
% 2nd column to names in the data and the 3rd column to the measured values
% in the dataset 
map = zeros(size(modelgenes,1),3)*NaN;
map = num2cell(map);
map(:,1) = modelgenes;
[~,idx] = ismember(modelgenes,table2cell(datanames));
map(:,2) = num2cell(idx);

% fill in corresponding transcription data values 
datavalues = table2cell(datavalues);
for i = 1:size(map,1)
    if ~isequal(map{i,2},0)
        map(i,3) = datavalues(map{i,2}); % FEHLER FIX THIS 
    end
end

% sets non measured at-codes to their original name (instead of the nan) -
% needed because later the at-codes will be replaced by the values based on
% the mapping, resulting in no exchange if the values are subsituted back
% to the original name 
for i = 1:size(map,1)
    if isequal(map{i,2},0)
        map{i,3} = map{i,1};
    end
end

%%%%%%%%%
%%% 4 %%%
%%%%%%%%%
% 4) replace genenames with their corresponding measured value (if not
% measured, the name genename stays in the rule)
for i=1:size(map,1)
    for j=1:length(rules)
        rules(j)=strrep(rules(j),map{i,1},num2str(map{i,3}));
    end
end
% find better solution, very slow  



%%%%%%%%%
%%% 5 %%%
%%%%%%%%%
% 5) store rules on seperate object, to save 'rules' 
mrm=rules;
modelrules=rules; %safety

%%%%%%%%%
%%% 6 %%%
%%%%%%%%%
% 6) Exclude non measured atcodes from the calculation - not entirely 
%    finished, might not account for every possible setting  
if isequal(excludenonmeasured,1)
    
    % find rules  
    mappingnumbers = map(:,3) ;
    strID = ~cellfun(@(x) isnumeric(x), mappingnumbers);
    unresolved_rules = mappingnumbers(strID);
    atID = cellfun(@(x) startsWith(x,'AT'),unresolved_rules);
    o = string(unresolved_rules(atID)); % o for simplicity, corresponds to unresolved ruleswithatgs
    %atrules = mrm(contains(mrm,['AT']));
    atrules = mrm;

    pattern = ["12*", "2*", "3*", "4*", "6*", "8*"];

    for rule = 1:length(atrules)
        for atg = 1:size(o,1)
          if contains(atrules(rule),o(atg))  
            if contains(atrules(rule),strcat('&  ',{' '},o(atg)))
                    atrules{rule} = char(strrep(atrules{rule},strcat('& ',{' '},o(atg)),''));

            elseif contains(atrules(rule),strcat('| ',{' '},o(atg)))
                    atrules{rule} = char(strrep(atrules{rule},strcat('| ',{' '},o(atg)),''));

            elseif contains(atrules(rule),strcat('(',o(atg),{' '},'&'))
                    atrules{rule} = char(strrep(atrules{rule},strcat('(',o(atg),{' '},'&'),'('));

            elseif contains(atrules(rule),strcat('(',o(atg),{' '},'|'))
                    atrules{rule} = char(strrep(atrules{rule},strcat('(',o(atg),{' '},'|'),'('));

            elseif contains(atrules(rule),strcat('(',{' '},o(atg),{' '},'|'))
                    atrules{rule} = char(strrep(atrules{rule},strcat('(',{' '},o(atg),{' '},'|'),'('));

            elseif contains(atrules(rule),strcat('(',{' '},o(atg),{' '},'&'))
                    atrules{rule} = char(strrep(atrules{rule},strcat('(',{' '},o(atg),{' '},'&'),'('));

            elseif contains(atrules(rule),'AT')
                    for patternID = 1:length(pattern)
                        if contains(atrules(rule),strcat(pattern(patternID),'(',o(atg),')',{' '},'&')) 
                            atrules{rule} = char(strrep(atrules{rule},strcat(pattern(patternID),'(',o(atg),')',{' '},'&'),''));
                        elseif contains(atrules(rule),strcat(pattern(patternID),'(',o(atg),')',{' '},'|'))
                            atrules{rule} = char(strrep(atrules{rule},strcat(pattern(patternID),'(',o(atg),')',{' '},'|'),''));
                        end
                    end
            else
                disp([atg,rule]);
            end
          end
        end
    end
    mrm = atrules; 
end



%%%%%%%%%
%%% 7 %%%
%%%%%%%%%
%%% 7) Solve rules 
% calculate values based on rules 
for n=1:length(mrm) 
    % Selects substring within last pair of brackets
    
while isempty(strfind(mrm{n},'('))==0

    lastbracket=strfind(mrm{n},'(');
    firstbracket=strfind(mrm{n},')');
    lastsubstr=mrm{n}(max(lastbracket):firstbracket(find(firstbracket>max(lastbracket),1)));

if ~contains(lastsubstr,'&')==0
    lastsubstr=strrep(lastsubstr,'&',',');
    lastsubstr=strrep(lastsubstr,'(','[');
    lastsubstr=strrep(lastsubstr,')',']');
    lastsubstr=strcat('min','(',lastsubstr,')');
    
elseif ~contains(lastsubstr,'|')==0
    if modellogic==0
        lastsubstr=strrep(lastsubstr,'|',',');
        lastsubstr=strrep(lastsubstr,'(','[');
        lastsubstr=strrep(lastsubstr,')',']');
        lastsubstr=strcat('max','(',lastsubstr,')');
    elseif modellogic==1
        lastsubstr=strrep(lastsubstr,'|','+');
    end
end

     % Evaluates expression within last brackets and substitutes value in logical rule

    try
       lastsubeval{n}=eval(lastsubstr); 
       mrm{n}=strrep(mrm{n},mrm{n}(max(lastbracket):firstbracket(find(firstbracket>max(lastbracket),1))),sprintf('%d',lastsubeval{n}));
    catch
        lastsubeval{n}=0;
        mrm{n}=strrep(mrm{n},'(','');mrm{n}=strrep(mrm{n},')','');
    end
     
end

     % Performs substitution when no brackets left 

if isempty(strfind(mrm{n},'('))==1 && (isempty(strfind(mrm{n},'&'))==0)&&(isempty(strfind(mrm{n},'|'))==1)
   mrm{n}=strrep(mrm{n},'&',',');
   mrm{n}=strcat('min','(','[',mrm{n},']',')');   
    
elseif isempty(strfind(mrm{n},'('))==1 &&(isempty(strfind(mrm{n},'|'))==0)&&(isempty(strfind(mrm{n},'&'))==1)
    if modellogic==0
        mrm{n}=strrep(mrm{n},'|',',');
        mrm{n}=strcat('max','(','[',mrm{n},']',')');
    elseif modellogic==1
        mrm{n}=strrep(mrm{n},'|','+');
    end
elseif isempty(strfind(mrm{n},'('))==1 && (isempty(strfind(mrm{n},'&'))==0)&&(isempty(strfind(mrm{n},'|'))==0)
    mrm{n}=strrep(mrm{n},'&',',');
    isor=strfind(mrm{n},'|');
    if length(isor)==1
    mrm{n}=strcat('min','(','[',mrm{n}(1:isor(1)-1),']',')',mrm{n}(isor(1)),'min','(','[',mrm{n}(isor(1)+1:length(mrm{n})),']',')');
    elseif length(isor)>1
        mrm{n}=strcat('min','(','[',mrm{n}(1:isor(1)-1),']',')',mrm{n}(isor(1):length(mrm{n})));
        for q=2:length(isor)
            if q==length(isor)
               isor=strfind(mrm{n},'|');
               mrm{n}=strcat(mrm{n}(1:isor(q-1)),'min','(','[',mrm{n}(isor(q-1)+1:isor(q)-1),']',')',mrm{n}(isor(q)),'min','(','[',mrm{n}(isor(q)+1:length(mrm{n})),']',')');
            else
              isor=strfind(mrm{n},'|');
              mrm{n}=strcat(mrm{n}(1:isor(q-1)),'min','(','[',mrm{n}(isor(q-1)+1:isor(q)-1),']',')',mrm{n}(isor(q):length(mrm{n})));
            end            
        end
        
    end

end
end
  

for n=1:length(mrm)
    if isempty(mrm{n})==0
        if modellogic==0
        mrm{n}=strrep(mrm{n},'|',',');
        mrm{n}=strcat('max','(','[',mrm{n},']',')');
        elseif modellogic==1
        mrm{n}=strrep(mrm{n},'|','+');
        end
        try
        mrm{n}=strrep(mrm{n},mrm{n},sprintf('%d',eval(mrm{n})));
        catch
            mrm{n}=modelrules{n};
        end
    end
end

% #*atcode - gets calculated!!! - so for example 12*atg123 gets calculated
calculated_rules=str2double(mrm);

%%%%%%%%%
%%% 8 %%%
%%%%%%%%%
% normalizes data to maximum value found if specified
if normalize==1
    calculated_rules = str2double(mrm);
    calculated_rules(:) = calculated_rules(:)/(max(calculated_rules));
    calculated_rules(isnan(calculated_rules))=1;
    vmax_eflux = calculated_rules;
end

%%%%%%%%%
%%% 9 %%%
%%%%%%%%%
% 9) Unresolved rules get set to their original ub/lb value 
if keepmodel==1
    vmaxeflux=calculated_rules;
    vmaxeflux(isnan(vmaxeflux))=vmax(isnan(vmaxeflux));
    calculated_rules = vmaxeflux;
end

%%%%%%%%%
%%% 10 %%
%%%%%%%%%
% 10) checks arbitrarily if everything was calculated correctly and sets
%     vmax to the calculated values 
if (isequal(size(calculated_rules),size(vmax)))
    vmax_eflux = calculated_rules;
else
    disp('mistake in sizes of rules/ub');
end



% comments:
%
% some at-codes are not in the dataset but in the rules
% exchange their values by the max measured value BUT BE CAREFULE BECAUSE 
% MULTIPLICATIVE NATURE OF THE PROBLEM
% OR JUST THINK OF A WAY ON HOW TO SOLVE THIS PROBLEM 
% maybe let everything that couldnt be calculated be 1 since the normation
% has to be done at some point anyways
%
%
% how does the normalization to 1 affect the overall calculations?
%
% use the changes in names done in local script here to main script for
% correct use and calculations 

end 



