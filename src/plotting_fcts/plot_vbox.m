function h = plot_vbox(x,values,color,quantiles,opt)
% h = plot_vbox(x,values,color,quantiles,opt)
%   quantiles; default = [.05 .25 .5 .75 .95];
%   color; default = Plotting_parameters.gray;
%   opt = {'linewidth' 'xwidth' 'Outcolor' 'BoxLine'}

global Plotting_parameters
Generate_Plotting_parameters

if isempty(values)
    warning('plot_vbox : empty value matrix; over')
    return
end

if isrow(values)
    values=values';
elseif all(size(values)>2) && ismatrix(values) && any(length(x)==size(values))
    ih = ishold;
    if ~exist('color','var')
        color = [];quantiles=[];opt=[];
    elseif ~exist('quantiles','var')
        quantiles=[];opt=[];
    elseif ~exist('opt','var')
        opt = [];
    end
    h = NaN(5,length(x));
    if length(x)~=size(values,2)
        values = values';
    end
    for i=1:length(x)
        h(:,i) = plot_vbox(x(i),values(:,i),color,quantiles,opt);
        hold on
    end
    if ~ih
        hold off
    end
    return
end

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
if ~exist('color','var') || isempty(color)
color = Plotting_parameters.gray;
end

Outcolor = Plotting_parameters.gray;
linewidth = 1.5;
xwidth = .4;
BoxLine = 'none';
if exist('opt','var')
    vars = {'linewidth' 'xwidth' 'Outcolor' 'BoxLine'};
    for i=1:length(vars)
        if isfield(opt,vars{i})
            eval([vars{i}  ' = opt.' vars{i} ';'])
        end
    end
end

ih = ishold;

Gq = quantile(values,quantiles);

h(3) = plot([1 1]*x, Gq([1 5]),'-','color',color,'linewidth',linewidth);
hold on
if Gq(4)~=Gq(2)
    h(2) = rectangle('position',[x-xwidth/2 Gq(2) xwidth Gq(4)-Gq(2)],'facecolor',color,'linestyle',BoxLine);
end
h(1) = line(x+[-.5 .5]*xwidth, Gq(3)*[1 1], 'color','k','linewidth',linewidth);
j = find( values<Gq(1) | values>Gq(5))';
if ~isempty(j)
    h(4) = plot(x*ones(size(j)),values(j),'.','color',Outcolor,'markersize',8);
end


h(5) = plot(NaN,NaN,'-','color',color,'linewidth',2*linewidth);

if ~ih
    hold off
end
