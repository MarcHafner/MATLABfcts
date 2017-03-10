function handle = addaxesletter(letter, axishandle, posshift, varargin)
% handle = addaxesletter(letter, axishandle)
%   add letter in upper left corner

if ~exist('axishandle','var') || isempty(axishandle)
    axishandle = gca;
else
    oldaxis = gca;
    set(gcf,'currentaxes', axishandle)
end

if ~exist('posshift','var') || isempty(posshift)
    posshift = [0 0];
end
pos = get(axishandle,'position');
h = annotation('textbox', ...
    [max(0,pos(1)-.05+posshift(1)) min(.95,pos(2)+pos(4)+posshift(2)) ...
    min(.05,pos(1)+posshift(1)) .05], ...
    'string', letter, 'linestyle', 'none', 'fontweight','bold', ...
    'fontsize', 12, 'horizontalalign', 'right', 'verticalalign', 'bottom', ...
    varargin{:});

if exist('axishandle','var') && ~isempty(axishandle)
    set(gcf,'currentaxes', oldaxis)
end

if nargout==1
    handle = h;
end