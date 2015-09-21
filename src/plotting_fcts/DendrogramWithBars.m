function [axis_handles, bar_handles, outperm, groups, dendo_handles] = ...
    DendrogramWithBars(SampleDist, t_SampleLabels, AxisPos, xwidth, ...
  PlotOptions, colorthreshold, linkagetype)
% function [axis_handles, bar_handles, outperm, groups, dendo_handles] =
%       DendrogramWithBars(SampleDist, t_SampleLabels, AxisPos, xwidth, ...
%           PlotOptions, colorthreshold, linkagetype)
%
%   if colorthreshold is negative, set as quantile of the distances.

labels = cellstr2str( table2cellstr(t_SampleLabels,0));

if ~isvector(SampleDist)
    SampleDist = squareform(SampleDist,'tovector');
end

if ~exist('linkagetype','var')
linkagetype = 'average';
end

axis_handles = get_newaxes(AxisPos);
tree = linkage(SampleDist,linkagetype);
if ~exist('colorthreshold','var') || isempty(colorthreshold)
    [dendo_handles,~,outperm] = dendrogram(tree,0,'labels',labels,'orientation','left');
    for i=1:length(dendo_handles), set(dendo_handles(i),'color','k'),end
else
    if colorthreshold<0
        colorthreshold = quantile(tree(:,3), -colorthreshold);
    end
    [dendo_handles,~,outperm] = dendrogram(tree,0,'labels',labels,'orientation','left',...
        'colorthreshold',colorthreshold);
    SampleColor = NaN(length(labels),3);
    for i=1:length(dendo_handles)
        idx = find(get(dendo_handles(i),'xdata')==0);
        if ~isempty(idx)
            ypos = get(dendo_handles(i),'ydata');
            for j=1:length(idx)
                if mod(ypos(idx(j)),1)==0
                    SampleColor(ypos(idx(j)),:) = get(dendo_handles(i),'color');
                end
            end
        end
    end
end

ylabs = get(gca,'yticklabel');
set(gca,'fontsize',8, 'ytick',[])
ylim([.5 length(labels)+.5])
xlim([-.2 1.1])

%%
fields = varnames(t_SampleLabels);
bar_handles = cell(length(fields),1);
groups = ones(height(t_SampleLabels),1);
cnt = 1;
for iF = 1:length(fields)
    axis_handles(iF+1) = get_newaxes([AxisPos(1)+AxisPos(3)+(iF-1)*xwidth AxisPos(2) xwidth AxisPos(4)],1);

    temp = t_SampleLabels.(fields{iF});
    temp = temp(outperm);

    for i=1:length(PlotOptions.(fields{iF}))
        idx = find(strcmp(cellstr(temp),PlotOptions.(fields{iF})(i)));
        bar_handles{iF}(i) = barh([-2;-1;idx], ones(size(idx,1)+2,1), 1, 'facecolor', ...
            PlotOptions.([fields{iF} 'Colors'])(i,:), 'linestyle','none');
    end

    if exist('SampleColor','var')
        for i = find(any(diff(SampleColor,1)~=0,2))'
            plot([0 1], [i i]+.5, '-k')
            cnt = cnt+1;
            groups((i+1):end) = cnt;
        end
    end

    ylim([.5 length(labels)+.5])
    xlim([.1 .9])
    set(gca,'xtick',[],'ytick',[])
end
set(gca,'ytick',1:length(labels), 'yticklabel', ylabs, 'fontsize',6, ...
    'yaxislocation','right','xtick',[])
