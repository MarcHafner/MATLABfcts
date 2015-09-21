function varargout = plot_violin(x, y, ybins, width, varargin)
% h = plot_violin(x, y, ybins, width, varargin)
%

if length(x)>1
    assert(any (length(x) == size(y)), 'x and y must have a common dimension')
    if length(x) ~= size(y,1)
        idx = find(length(x) == size(y),1,'first');
        y = permute(y,[ idx setdiff(1:length(size(y)),idx)]);
    end

    if ~exist('ybins','var') || isempty(ybins)
        ydiff = (max(y(:))-min(y(:)))/20;
        ybins = (min(y(:))-ydiff*3):( ydiff ):(max(y(:))+ydiff*3);
    end

    if ~exist('width','var') || isempty(width)
        width = .7;
    end

    ish = ishold;
    hold on
    for i=1:length(x)
        h = plot_violin(x(i), y(i,:), ybins, width, varargin{:});
    end
    if ~ish
        hold off
    end
    return
end

if ~exist('width','var') || isempty(width)
    width = .7;
end

if ~exist('ybins','var') || isempty(ybins)
    ydiff = (max(y)-min(y))/20;
    ybins = (min(y)-ydiff*3):( ydiff ):(max(y)+ydiff*3);
else
    ydiff = diff(ybins(1:2));
end

ydist = ksdensity(y, ybins, 'width', ydiff*1.3);
range = max(find(ydist>1e-6*max(ydist),1,'first')-2,1):...
    min(find(ydist>1e-6*max(ydist),1,'last')+2,length(ydist));

ydist = ToRow(ydist)*width/2/max(ydist);


ish = ishold;
hold on

% h = plot([x+ydist(range) NaN x-ydist(range)], [ybins(range) NaN ybins(range)], ...
%     varargin{:});

h = patch([x+ydist(range) x x-ydist(range(end:-1:1)) x], ...
    ybins(range([1:end end end:-1:1 1])), 'k', ...
    varargin{:});

h(2) = plot(x, nanmedian(y), '.k');

if nargout>0
    varargout = {h};
end

if ~ish
    hold off
end
