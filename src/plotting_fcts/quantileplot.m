function h = quantileplot(T, data, color, quantiles, alpha)
% h = quantileplot(T, data, color, quantiles, alpha)
%   plot the quantiles of a trajectory with transparancy alpha

ish = ishold;

if ~exist('quantiles','var') || isempty(quantiles)
    quantiles = [.1 .25 .5  .75 .9];
end
if ~exist('color','var') || isempty(quantiles)
    color = [.6 .6 .6];
end

if ~exist('alpha','var')
    alpha = .7;
end

interquart = quantile(data,quantiles([2 4]),2)-[zeros(size(data,1),1) quantile(data,quantiles(2),2)];
first = find(any(~isnan(interquart),2), 1, 'first');
last = find(any(~isnan(interquart),2), 1, 'last');
h = area(T(first:last), interquart(first:last,:) );
set(h(1),'facecolor','none')
set(h(2),'facecolor',color,'edgecolor','k')
set(get(h(2),'children'),'facealpha',alpha)
hold on
h(3:4) = plot(T, quantile(data,quantiles([1 5]),2), '--','color',color,'linewidth',1);
h(5) = plot(T, quantile(data,quantiles(3),2), '-k','linewidth',2);

if ~ish
    hold off
end
