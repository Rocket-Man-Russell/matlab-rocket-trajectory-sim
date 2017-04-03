function [varargout] = extendSimData(expFlightDur,timeStep,simRunNum,numOfSteps,varargin)
%
if      exist('wextend.m','file') == 0
    
        for count = 3:nargin
        [varargout{count-2}] = deal(zeros(1,(expFlightDur/timeStep)));
        end
    
elseif  exist('wextend.m','file') == 2
    
    if simRunNum > 1
    
        for count = 3:nargin
        [varargout{count-2}] = 0*wextend('addcol','zpd',[varargin{count-2}],(numOfSteps-length(varargin{count-2})),'r');
        end
    end
end
%