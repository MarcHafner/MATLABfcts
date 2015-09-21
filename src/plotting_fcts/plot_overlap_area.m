function [ha, hl] =plot_overlap_area(x, y1, y2, colors, varargin)
% plot_overlap_area(x, y1, y2, colors, varargin)
% y1 and y2 are Nx2 matrices with upper and lower bound of the area
% y1 and y2 variables should have the same x sampling (N points).
% assume that y(:,1) <= y(:,2) for y1 and y2
%
% may not handle well multiple overlapping area(for now)
%

if ~exist('colors','var') || isempty(colors)
    colors = [1 .5 .4; .4 .5 1; .6 .6 .6];
end
linecolors = max(colors.^.8-.1,0);

% define the intersection
y_inter = NaN(length(x),2);
for i=1:length(x)
    if ~(all(y2(i,1)>y1(i,:)) || all(y2(i,2)<y1(i,:)) || ...
            all(y2(i,:)>y1(i,2)) || all(y2(i,:)<y1(i,1)))
        y_inter(i,:) = [max(y1(i,1),y2(i,1)) min(y1(i,2),y2(i,2))];
    end
end

idx = find(~isnan(y1(:,1)));
p1 = [ToColumn(x([idx(1:end); idx(end:-1:1); idx(1)])) ...
    [y1(idx(1:end),1); y1(idx(end:-1:1),2); y1(idx(1),1)]];

idx = find(~isnan(y2(:,1)));
p2 = [ToColumn(x([idx(1:end); idx(end:-1:1); idx(1)])) ...
    [y2(idx(1:end),1); y2(idx(end:-1:1),2); y2(idx(1),1)]];

idx = find(~isnan(y_inter(:,1)));
pi = [ToColumn(x([idx(1:end); idx(end:-1:1); idx(1)])) ...
    [y_inter(idx(1:end),1); y_inter(idx(end:-1:1),2); y_inter(idx(1),1)]];


ish = ishold;
hold on

ha = patch(p1(:,1), p1(:,2), colors(1,:));
ha(2) = patch(p2(:,1), p2(:,2), colors(2,:));
ha(3) = patch(pi(:,1), pi(:,2), colors(3,:));

hl = line(p1(:,1), p1(:,2), 'color', linecolors(1,:), varargin{:});
hl(2) = line(p2(:,1), p2(:,2), 'color', linecolors(2,:), varargin{:});
hl(3) = line(pi(:,1), pi(:,2), 'color', linecolors(3,:), varargin{:});
