function makeGraph(type,simRunNum,graphTitle,xLabel,yLabel,varargin)

plotVars=[];
legendVars=[];
figure


for count = 1:length(varargin)
    if      ~(iscellstr(varargin(count)))
            plotVars = [plotVars;varargin(count)];
            
    elseif  iscellstr(varargin(count))
            legendVars = [legendVars;varargin(count)];
            
    else
            error('%s\n','Input argument is not a valid type for graph production.', ...
                'Must be either character array for legends or numeric array for x/y point ranges.')
    end
end


[plotVarRows,plotVarCols] = size(plotVars);
if (rem(plotVarRows,1) == 0)
    
    for count = 1:2:plotVarRows
        if      strcmp(type,'Line') == 1
                plot(plotVars{count},plotVars{count+1})
        elseif  strcmp(type,'Scatter') == 1
                scatter(plotVars{count},plotVars{count+1})
        end
        hold on
    end
    
elseif ~(rem(plotVarRows,1) == 0)
        error('Graph plotting requires pairs of vectors for data to plot, but an odd number have been provided.')
end


title(graphTitle)
xlabel(xLabel)
ylabel(yLabel)
xlim([0,inf])
   

if ~(isempty(legendVars))
    [legendVarRows,legendVarCols] = size(legendVars);
    
    if ((plotVarRows/legendVarRows) == 2)
        legend(legendVars)
    else
        warning('An incorrect number of legends have been entered for the number of graphs being plotted.')
    end
    
end


saveas(gcf,[graphTitle,'_(',num2str(simRunNum),').jpg'])
graphDelay = 1.5;
tic; while toc < graphDelay; end
close
