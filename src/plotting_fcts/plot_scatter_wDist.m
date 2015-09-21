function [a,h] = plot_scatter_wDist(xvals, yvals, colors, outerpos, varargin)
% [a,h] = plot_scatter_wDist(xvals, yvals, colors, outerpos, varargin)
%
global Plotting_parameters

% process the data
if iscell(xvals)
    xvals = cellfun2(@ToRow, xvals);
elseif isvector(xvals)
    xvals = {ToRow(xvals)};
elseif ismatrix(xvals)
    temp = xvals;
    xvals = cell(1,size(xvals,2));
    for i=1:size(xvals,2)
        xvals{i} = temp(:,i)';
    end
end

if iscell(yvals)
    yvals = cellfun2(@ToRow, yvals);
elseif isvector(yvals)
    yvals = {ToRow(yvals)};
elseif ismatrix(yvals)
    temp = yvals;
    yvals = cell(1,size(yvals,2));
    for i=1:size(yvals,2)
        yvals{i} = temp(:,i)';
    end
end
assert(length(xvals)==length(yvals))
assert(all(cellfun(@(x,y) length(x)==length(y), xvals, yvals)))
nData = length(xvals);

if isstr(colors)
    if length(colors)==1
        colors = repmat(colors, nData,1);
    else
        assert(length(colors)==nData)
        colors = ToColumn(colors);
    end
else % RGB vecotr/matrix
    assert(size(colors,2)==3)
    if size(colors,1)==1
        colors = repmat(colors, nData,1);
    else
        assert(size(colors,1)==nData)
    end
end

if ~exist('outerpos','var')
    outerpos = [.1 .1 .85 .85];
end
assert(all(outerpos<1))

% get the varargin and default
p = inputParser;
addParameter(p,'xplot_height',outerpos(4)*.125, @(x) x<outerpos(4)*.5);
addParameter(p,'yplot_width',outerpos(3)*.125, @(x) x<outerpos(3)*.5);
addParameter(p,'xcenters',[], @isvector);
addParameter(p,'ycenters',[], @isvector);
addParameter(p,'xwidth',[], @isvector);
addParameter(p,'ywidth',[], @isvector);
addParameter(p,'markers','o',@(x) ischar(x) && ismember(length(x),[1 nData]));
parse(p,varargin{:});
p = p.Results;

if length(p.markers)==1, p.markers=repmat(p.markers,1,nData); end

if isempty(p.xcenters)
    p.xcenters = [quantile([xvals{:}],.01) quantile([xvals{:}],.99)] ;
    p.xcenters = p.xcenters(1):(diff(p.xcenters)/200):p.xcenters(2);
end
if isempty(p.xwidth)
    p.xwidth = 1.5*diff(p.xcenters(1:2));
end

if isempty(p.ycenters)
    p.ycenters = [quantile([yvals{:}],.01) quantile([yvals{:}],.99)] ;
    p.ycenters = p.ycenters(1):(diff(p.ycenters)/200):p.ycenters(2);
end
if isempty(p.ywidth)
    p.ywidth = 1.5*diff(p.ycenters(1:2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% scatter plot
innerpos = outerpos+[p.yplot_width*1.1 p.xplot_height*1.1 ...
    -p.yplot_width*1.1 -p.xplot_height*1.1];
a = get_newaxes(innerpos,1);
for i=1:nData
    h(i) = plot(xvals{i}, yvals{i}, p.markers(i), 'color', colors(i,:));
end
xlims = xlim;
ylims = ylim;
set(gca, Plotting_parameters.axes{:}, 'xtick', [], 'ytick', [])

% x-axis dist
a(2) = get_newaxes([outerpos(1:2)+[p.yplot_width*1.1 0] innerpos(3) p.xplot_height],1);

maxf = 0;
for i=1:nData
    f = ksdensity(xvals{i}, p.xcenters, 'width', p.xwidth);
    plot(p.xcenters, f, 'color', colors(i,:));
    maxf = max([f';maxf]);
end
xlim(xlims)
ylim([0 1.1*maxf])
set(gca, Plotting_parameters.axes{:}, 'ytick', [])

% x-axis dist
a(3) = get_newaxes([outerpos(1:2)+[0 p.xplot_height*1.1] p.yplot_width innerpos(4)],1);
maxf = 0;
for i=1:nData
    f = ksdensity(yvals{i}, p.ycenters, 'width', p.ywidth);
    plot(f, p.ycenters, 'color', colors(i,:));
    maxf = max([f';maxf]);
end
xlim([0 1.1*maxf])
ylim(ylims)
set(gca, Plotting_parameters.axes{:}, 'xtick', [])
