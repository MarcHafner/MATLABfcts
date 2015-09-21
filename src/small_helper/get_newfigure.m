function fhandle = get_newfigure(number, position, varargin)
% fhandle = get_newfigure(number, position, varargin)
%   get a clean new figure
%   first argument of varargin can be the filename
%   further should be pairs of parameter/value as in set(gcf, ...)
%

if ~exist('position','var') || isempty(position)
    position = [100 80 600 450];
end

if nargout==0
    figure(number)
else
    fhandle = figure(number);
end
clf

if (length(varargin)==1 || mod(length(varargin),2)==1) && varargin{1}(end-3)=='.'
    varargin = ['filename' varargin];
end


set(gcf, 'color','w','paperunits','inches', 'PaperPositionMode', 'auto');
set(gcf,'papersize', [position(3)/90 position(4)/90]);
set(gcf, 'position', position, varargin{:});
