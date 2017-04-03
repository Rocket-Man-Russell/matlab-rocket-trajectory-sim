function [varName] = changeInput(simRun,varDescriptor,inputQuestion,varName,varargin)
%
global yes
%
if      simRun > 1
    
        if (max(strcmp((input(['Change input value for ',varDescriptor,'? [Y/N]: '],'s')),yes))) == true
            changeValue = true;
        else
            changeValue = false;
        end
        
elseif  simRun == 1
        changeValue = true;
end
%
if      changeValue == true
            validInput = false;

	while 	validInput == 0
		varName = input(inputQuestion,'s');
        varName = str2num(varName);

		if 	( (nargin > 4) && (isempty(find(varName == varargin{1})) == false) ) ||...
                ...
                (  (nargin == 4) && ...
                ~(isempty(varName)) && ...
			     isnumeric(varName) && ...
			     isreal(varName) && ...
			     isscalar(varName) == true  )
             
            validInput = true;
        else
            warning('You did not enter a valid input. Try again.')
		end
	end
        
elseif  changeValue == false
        varName = varName;
end
