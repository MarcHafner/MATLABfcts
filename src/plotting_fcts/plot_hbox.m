function h = plot_hbox(y,values,color,quantiles,opt)
% h = plot_hbox(y,values,color,quantiles,opt)
%   quantiles; default = [.05 .25 .5 .75 .95];
%   color; default = Plotting_parameters.gray;
%   opt = {'linewidth' 'xwidth' 'Outcolor' 'BoxLine'}

global Plotting_parameters
Generate_Plotting_parameters


if exist('color','var') && isstruct(color)
    opt = color;
    color = [];
elseif exist('quantiles','var') && isstruct(quantiles)
    opt = quantiles;
    quantiles = [];
end

if ~exist('quantiles','var') || isempty(quantiles)
    quantiles = [.05 .25 .5 .75 .95];
end
if ~exist('color','var')
color = Plotting_parameters.gray;
end

Outcolor = Plotting_parameters.gray;
linewidth = 1.5;
xwidth = .4;
if exist('opt','var')
    vars = {'linewidth' 'color' 'xwidth' 'Outcolor'};
    for i=1:length(vars)
        if isfield(opt,vars{i})
            eval([vars{i}  ' = opt.' vars{i}])
        end
    end
end

ih = ishold;

Gq = quantile(values,quantiles);

h(3) = plot(Gq([1 5]),[1 1]*y, '-','color',color,'linewidth',linewidth);
hold on
if Gq(4)~=Gq(2)
    h(2) = rectangle('position',[Gq(2) y-xwidth/2 Gq(4)-Gq(2) xwidth],'facecolor',color,'linestyle','none');
end
h(1) = line(Gq(3)*[1 1], y+[-.5 .5]*xwidth, 'color','k','linewidth',linewidth);
for j = find( values<Gq(1) | values>Gq(5))'
    temph = plot(values(j),y,'.','color',Outcolor);
    h(4) = temph(1);
end

if ~ih
    hold off
end
