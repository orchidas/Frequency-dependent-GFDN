function [plotAxes, plotHandles] = plotImpulseResponseMatrix( t, ir, varargin )
%plotImpulseResponseMatrix - Plot matrix of impulse response / polynomial matrix in subplots
%
% Syntax:  [plotAxes, plotHandles] = plotImpulseResponseMatrix( t, ir, varargin )
%
% Inputs:
%    t - x-values, e.g., time-domain variable of size [FIR,1] or [1,FIR]
%    ir - matrix of y-values of size [out, in, FIR]
%    varargin - extra plotting parameters
%
% Outputs:
%    plotAxes - Matrix of axes handles; size [out, in]
%    plotHandles - Matrix of axes handles; size [out, in]
%
% Example: 
%    [plotAxes, plotHandles] = plotImpulseResponseMatrix( 1:100, randn(3,2,100), 'xlabel','Time (seconds)','ylim',[-1, 1])
%
% Author: Dr.-Ing. Sebastian Jiro Schlecht, 
% Aalto University, Finland
% email address: sebastian.schlecht@aalto.fi
% Website: sebastianjiroschlecht.com
% 29. December 2019; Last revision: 5. April 2021
% Modified by: Orchisama Das, January 2022.

%% parse input
p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'xlabel',[]);
addParameter(p,'ylabel',[]);
addParameter(p,'title',[]);
addParameter(p,'xlim',[]);
addParameter(p,'ylim',[]);
addParameter(p, 'stemFlag',[]);
addParameter(p, 'colors', []);
addParameter(p, 'xlogaxis', 0);
addParameter(p,'save',0);
addParameter(p,'marker',[]);
addParameter(p, 'plotYLabels',1);
addParameter(p, 'commonYAxis',0);
parse(p,varargin{:});

xLabel = p.Results.xlabel;
yLabel = p.Results.ylabel;
Title = p.Results.title;
xLim = p.Results.xlim;
yLim = p.Results.ylim;
stemFlag = p.Results.stemFlag;
saveFlag = p.Results.save;
cols = p.Results.colors;
xlogaxis = p.Results.xlogaxis;
marker = p.Results.marker;
plotYLabels = p.Results.plotYLabels;
commonYAxis = p.Results.commonYAxis;
plotArg = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)];


%% plot
if isempty(t)
    t = 1:size(ir,3);
end

if isempty(stemFlag)
    stemFlag = 1;
end

% change order to [FIR, out, in]
ir = permute(ir,[3 1 2]);
numberOfOutputs = size(ir,2);
numberOfInputs = size(ir,3);
plotAxes = [];

for itOut = 1:numberOfOutputs
    for itIn = 1:numberOfInputs
        plotAxes(itOut,itIn) = subplot(numberOfOutputs, numberOfInputs, sub2ind([numberOfInputs,numberOfOutputs], itIn, itOut));
        hold on; grid on;
        % plotHandles(itOut,itIn) = reduce_plot(t, ir(:,itOut,itIn),varargin{:});
        if stemFlag
            if isempty(cols)
                plotHandles(itOut,itIn) = stem(t, ir(:,itOut,itIn),plotArg{:}, 'MarkerSize',2); 
            else
                plotHandles(itOut,itIn) = stem(t, ir(:,itOut,itIn),plotArg{:}, 'Color', cols, 'MarkerSize',2); 
            end
        else
            if isempty(cols)
                if isempty(marker)
                    if xlogaxis
                        plotHandles(itOut,itIn) = semilogx(t, ir(:,itOut,itIn),plotArg{:}, 'LineWidth', 0.7); 
                    else
                        plotHandles(itOut,itIn) = plot(t, ir(:,itOut,itIn),plotArg{:}, 'LineWidth', 0.7); 
                    end
                else
                    if xlogaxis
                        plotHandles(itOut,itIn) = semilogx(t, ir(:,itOut,itIn),plotArg{:},'Marker',marker,'MarkerSize',3); 
                    else
                        plotHandles(itOut,itIn) = plot(t, ir(:,itOut,itIn),plotArg{:}, 'Marker',marker,'MarkerSize',3); 
 
                    end
                end
            else
                if isempty(marker)
                    if xlogaxis
                        plotHandles(itOut,itIn) = semilogx(t, ir(:,itOut,itIn),plotArg{:}, 'Color', cols); 
                    else
                        plotHandles(itOut,itIn) = plot(t, ir(:,itOut,itIn),plotArg{:}, 'Color', cols); 
                    end
                else
                    if xlogaxis
                        plotHandles(itOut,itIn) = semilogx(t, ir(:,itOut,itIn),plotArg{:}, 'Marker', marker,'MarkerSize', 3, 'Color', cols); 
                    else
                        plotHandles(itOut,itIn) = plot(t, ir(:,itOut,itIn),plotArg{:}, 'Marker', marker, 'MarkerSize', 3, 'Color', cols); 
                    end
                end

            end

        end
    end    
end


%% remote ticks of inner subplots
if ~plotYLabels
    set(plotAxes(:,2:end), 'YTickLabel', []);
end
set(plotAxes(1:end-1,:), 'XTickLabel', []);



%% set limits
if isempty(xLim) 
    xLims = cell2mat(get(plotAxes,'XLim'));
    commonXLim = [min(xLims(:,1)),max(xLims(:,2))];
    set(plotAxes,'XLim',commonXLim);
else
    set(plotAxes,'XLim',xLim);
end 

if isempty(yLim) 
    if commonYAxis
        yLims = cell2mat(get(plotAxes,'YLim'));
        commonYLim = [min(yLims(:,1)),max(yLims(:,2))];
        set(plotAxes,'YLim',commonYLim);
    end
else
    set(plotAxes,'YLim',yLim);
end



%% Give common xlabel, ylabel and title to your figure
commonAx=axes(gcf,'visible','off');
commonAx.Title.Visible='on';
commonAx.XLabel.Visible='on';
commonAx.YLabel.Visible='on';
ylabel(commonAx,yLabel);
xlabel(commonAx,xLabel);
title(commonAx,Title);

if saveFlag && plotYLabels
    xticklabel = get(plotAxes, 'XTickLabel');
    yticklabel = get(plotAxes, 'YTickLabel');
    k = 1;
    for i = 1:size(plotAxes,1)
        set(plotAxes(end,i), 'XTickLabel', xticklabel{i*numberOfInputs}, 'FontSize', 6);
        for j = 1:size(plotAxes,2)
            set(plotAxes(j,i), 'YTickLabel', yticklabel{k}, 'FontSize', 6);
            k = k+1;
        end
    end
    set(commonAx, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');
end


%% ensure x-axis is logarithmic if chosen to be so
if xlogaxis
    set(plotAxes, 'xscale', 'log');
end




axis tight;




