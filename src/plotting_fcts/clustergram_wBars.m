function h = clustergram_wBars(data, rowlabels,  columnlabels, varargin)
% h = clustergram_wBars(data, rowlabels,  columnlabels, varargin)
%
%   options:
%       RowPDist
%       ColumnPDist
%       Linkage
%       rowbarwidth
%       colbarwidth
%       PO
%       outerpos
%       rowannotations
%       columnannotations
%       cmap
%



Generate_Plotting_parameters;

p = inputParser;
addParameter(p,'RowPDist','correlation')
addParameter(p,'ColumnPDist','correlation')
addParameter(p,'Linkage','Average')
addParameter(p,'rowbarwidth', .04)
addParameter(p,'colbarwidth', .05)
addParameter(p,'PO', Plotting_parameters)
addParameter(p,'outerpos', NaN)
addParameter(p,'rowannotations', table)
addParameter(p,'columnannotations', table)
addParameter(p,'cmap', [.3 .3 .2;Plotting_parameters.cmapRW])
addParameter(p,'clim', [])

parse(p,varargin{:});
p = p.Results;
PO = p.PO;

if ~strcmp(p.RowPDist,'none')
    RowDist = pdist(data, p.RowPDist);
    TreeRows = linkage(RowDist, p.Linkage);
end

if ~strcmp(p.ColumnPDist,'none')
    ColumnDist = pdist(data', p.ColumnPDist);
    TreeColumns = linkage(ColumnDist, p.Linkage);
end

rowannotations = TableToCategorical(p.rowannotations);
columnannotations = TableToCategorical(p.columnannotations);

Nrowannot = size(rowannotations,2);
Ncolannot = size(columnannotations,2);

%% plotting positions

if any(isnan(p.outerpos)) || ~all(p.outerpos<1)
    outerpos = [.015+p.rowbarwidth .16 .88 .83-p.colbarwidth];
else
    outerpos = p.outerpos;
end
innerpos = outerpos+[p.rowbarwidth*Nrowannot 0 ...
    -p.rowbarwidth*Nrowannot -p.colbarwidth*Ncolannot];



%% rows Dendogram
if ~strcmp(p.RowPDist,'none')
    get_newaxes([outerpos(1)-p.rowbarwidth-.015 innerpos(2) p.rowbarwidth+.01 innerpos(4)])
    [h_RowDend,~,permRows] = dendrogram(TreeRows,0,'orientation','left');
    set(gca,'xtick',[],'ytick',[],'fontweight','bold','fontsize',8,'visible','off')
    ylim([.5 length(permRows)+.5])
else
    h_RowDend = [];
    permRows = 1:size(data,1);
end

% Row color bars
h_RowBars = cell(size(rowannotations,2),1);
for iR = 1:size(rowannotations,2)
    get_newaxes([outerpos(1)+(iR-1)*p.rowbarwidth innerpos(2) p.rowbarwidth-.005 innerpos(4)],1)

    h_RowBars{iR} = [];
    cases = unique(rowannotations.(iR));

    if isfield(PO, rowannotations.Properties.VariableNames{iR}) && ...
            isfield(PO, [rowannotations.Properties.VariableNames{iR} 'Colors'])
        colors = NaN(length(cases), 3);
        for i = 1:length(cases)
            if any(PO.(rowannotations.Properties.VariableNames{iR})==cases(i))
                colors(i,:) = PO.([rowannotations.Properties.VariableNames{iR} 'Colors'])(...
                    PO.(rowannotations.Properties.VariableNames{iR})==cases(i),:);
            else
                colors(i,:) = .95*[1 1 1];
            end
        end
    else
        colors = gray(length(cases));
    end

    for i = 1:length(cases)
        idx = find(rowannotations.(iR)(permRows)==cases(i));
        h_RowBars{iR}(i) = barh([-2;-3;idx], ones(length(idx)+2,1), 1, ...
            'linestyle','none','facecolor', colors(i,:));
    end

    set(gca,'ytick',[],'xtick',.5,'xticklabel',rowannotations.Properties.VariableNames{iR}, ...
        PO.axislabel{:},'xticklabelrotation',90)
    xlim([.2 .8])
    ylim([.5 length(permRows)+.501])
end



%% column dendrogram


if ~strcmp(p.ColumnPDist,'none')
    get_newaxes([innerpos(1) outerpos(2)+outerpos(4) innerpos(3) p.colbarwidth+.01])
    [h_ColDend,~,permCols] = dendrogram(TreeColumns,0,'orientation','top');
    set(gca,'xtick',[],'ytick',[],'fontweight','bold','fontsize',8,'visible','off')
    xlim([.5 length(permCols)+.5])
else
    h_ColDend = [];
    permCols = 1:size(data,2);
end

% Column color bars
h_ColBars = cell(size(columnannotations,2),1);
for iR = 1:size(columnannotations,2)
    get_newaxes([innerpos(1) outerpos(2)+outerpos(4)-iR*p.colbarwidth+.005 innerpos(3) p.colbarwidth-.005],1)

    h_ColBars{iR} = [];
    cases = unique(columnannotations.(iR));

    if isfield(PO, columnannotations.Properties.VariableNames{iR}) && ...
            isfield(PO, [columnannotations.Properties.VariableNames{iR} 'Colors'])
        colors = NaN(length(cases), 3);
        for i = 1:length(cases)
            if any(PO.(columnannotations.Properties.VariableNames{iR})==cases(i))
                colors(i,:) = PO.([columnannotations.Properties.VariableNames{iR} 'Colors'])(...
                    PO.(columnannotations.Properties.VariableNames{iR})==cases(i),:);
            else
                colors(i,:) = .95*[1 1 1];
            end
        end
    else
        colors = gray(length(cases));
    end

    for i = 1:length(cases)
        idx = find(columnannotations.(iR)(permCols)==cases(i));
        h_ColBars{iR}(i) = bar([-2;-3;idx], ones(length(idx)+2,1), 1, ...
            'linestyle','none','facecolor', colors(i,:));
    end

    set(gca,'xtick',[],'ytick',.5,'yticklabel',columnannotations.Properties.VariableNames{iR}, ...
        PO.axislabel{:})
    ylim([.2 .8])
    xlim([.5 length(permCols)+.501])
end


%% main heatmap
h_hmapax = get_newaxes(innerpos);

if isempty(p.clim)
    clim = [-1 1]*quantile(abs(data(~isnan(data))),.95);
else
    clim = p.clim;
end
temp = data;

if any(isnan(temp(:)))
    temp(~isnan(temp)) = max(min(clim)+diff(clim)/size(p.cmap,1), temp(~isnan(temp)));
end

h_map = imagesc(temp(permRows,permCols), clim);
colormap(p.cmap)


set(gca,PO.axes{:},'box','on','ydir','normal','yaxislocation','right')
if ~isempty(columnlabels)
    set(gca,'xtick',1:length(permCols),'xticklabel',columnlabels(permCols), 'xticklabelrotation',90)
else
    set(gca,'xtick',[])
end
if ~isempty(rowlabels)
    set(gca, 'ytick',1:length(permRows),'yticklabel',rowlabels(permRows))
else
    set(gca, 'ytick',[])
end


xlim([.5 length(permCols)+.5])
ylim([.5 length(permRows)+.501])

if nargout>0
    h = struct;
    h.RowDend = h_RowDend;
    h.RowBars = h_RowBars;
    h.ColDend = h_ColDend;
    h.ColBars = h_ColBars;
    h.heatmap = h_map;
    h.hmapax  = h_hmapax ;
end
