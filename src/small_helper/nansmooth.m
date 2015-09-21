function smoothed = nansmooth(input,step,method)
% smoothed = nansmooth(input,step,method)
%   rewritten smooth function to discard NaN

if ~exist('method','var')
    method = 'rlowess';
end
if ~exist('step','var')
    step = 5;
end

smoothed = NaN(size(input));
for j=1:size(input,2)
    first = find(~isnan(input(:,j)),1,'first');
    if isempty(first), first=1;end
    last = find(isnan(input(first:end,j)),1,'first')+first-2;
    if isempty(last), last=size(input,1);end
    smoothed(first:last,j) = smooth(input(first:last,j),step,method);
end
