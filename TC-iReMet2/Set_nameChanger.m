function [oi] = Set_nameChanger(oi,rxnpathway)


	%%% 
	%%% dummy is ranking.withTnA_man for example
	%%%
	
	dummy = cell2mat(oi(:,1:2));
	newnames = cell(size(dummy,1),1);
	for row = 1:size(dummy,1)
	
		% set checkpattern
		endingpattern = ["m","c","p","h"];
		endstr = rxnpathway(dummy(row,2),2);
		endstr = extractAfter(endstr,'_');
		% exchange with newname if according to compartment
		if any(strcmp(endingpattern,endstr))	
			comp = endingpattern(strcmp(endingpattern,endstr));
			bigname = rxnpathway{dummy(row,2),3};
			newname = strcat(bigname,{''},'[',comp,']');		
		else 
			newname = string(rxnpathway{dummy(row,2),3});
		end	
		% set new name 
		newnames{row} = newname;	
	end
	
	oi(:,3) = newnames;


end
















