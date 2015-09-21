
function a = categories_bars(t_group_data, fields, plotopt, ClusterStats, maxNumcat)
% a = categories_bars(t_group_data, fields, plotopt, ClusterStats)
%
%

Groupfiels = fields{1};
Gidx = t_group_data.(Groupfiels);
fields = fields(2:end);

if ~isfield(plotopt, 'labels')
    plotopt.labels = {'fontsize', 8, 'fontweight', 'bold'};
end
if ~isfield(plotopt, 'axis')
    plotopt.axis = {'fontsize', 6};
end

if exist('ClusterStats','var') && ~isempty(ClusterStats)
    Nclust = ClusterStats.Nclust;
    Cluster = 1;
    Nbars = length(ClusterStats.cnt);
else
    Nclust = max(Gidx);
    Cluster = 0;
    Nbars = max(Gidx);
end

if ~exist('maxNumcat','var')
    maxNumcat = 10;
end

%% distributions

Fiscat = false(length(fields),1);  % field is categorical
Fcats = cell(length(fields),1);  % field categories
Fsize = NaN(length(fields),1); % number of categories
Ffraction = cell(length(fields),1);
for iF = 1:length(fields)
    if iscategorical(t_group_data.(fields{iF}))
        Fiscat(iF) = true;
        if any(isundefined(t_group_data.(fields{iF})))
            warning('some undefined values for field %s', fields{iF})
            t_group_data.(fields{iF})(isundefined(t_group_data.(fields{iF}))) = ...
                'undef';
        end
        Fcats{iF} = unique(t_group_data.(fields{iF}));
        Fsize(iF) = length(Fcats{iF});
    elseif isnumeric(t_group_data.(fields{iF})) && ...
            length(unique(t_group_data.(fields{iF})))<=maxNumcat
        % maxNumcat or less values is considered category
        Fiscat(iF) = true;
        Fcats{iF} = unique(t_group_data.(fields{iF}));
        Fsize(iF) = length(Fcats{iF});
    elseif ~isnumeric(t_group_data.(fields{iF}))
        error('fields %s should be categorical or numeric', fields{iF})
    end

    if Fiscat(iF)
        isnum = isnumeric(Fcats{iF});
        if isnum
            Fcats{iF} = num2cellstr(Fcats{iF});
        else
            Fcats{iF} = cellstr(Fcats{iF});
        end
        if isfield(plotopt, [fields{iF} 'Colors']) && isfield(plotopt, fields{iF})
            Fcats{iF} = union(intersect(plotopt.(fields{iF}), Fcats{iF}, 'stable'), ...
                setdiff(Fcats{iF}, plotopt.(fields{iF}), 'stable'), 'stable');
        end

        Ffraction{iF} = zeros(Nbars, Fsize(iF));
        for i=1:Nbars
            if isnum
                Ffraction{iF}(i,:) = hist(t_group_data.(fields{iF})(Gidx==i),...
                    cellfun(@str2num,Fcats{iF}));
            else
                Ffraction{iF}(i,:) = hist_str(t_group_data.(fields{iF})(Gidx==i),...
                    Fcats{iF});
            end
        end
        Ffraction{iF} = NormSumUnit(Ffraction{iF},2);
    end
end
Fsize(isnan(Fsize)) = nanmean(Fsize); % for assigning the size of the plots



%%

xpos = .07;
xwidth = .72;
xlpos = .8;
barw = .7;
yspace = .02;

ypos = NaN(Cluster+length(fields)+1,2);
if Cluster
    ypos(1,:) = [.96  .035];
    ytot = .96;
else
    ytot = 1;
end
yunit = (ytot-(length(fields)+2.5)*yspace)/sum(Fsize);
ypos(end,:) = [yspace yspace*1.5];
for iF=length(fields):-1:1
    ypos(iF+Cluster,:) = [sum(ypos(iF+1+Cluster,:))+yspace Fsize(iF)*yunit];
end


get_newfigure(plotopt.fignumber, [50 50 650 750], [plotopt.figurename '.pdf']);

txtoutput = cell(2,1+Nbars);
txtoutput(1,2:(Nclust+1)) = num2cellstr(1:Nclust);
if Nclust==Nbars-1
    txtoutput(1,end) = {'Extra clust'};
end
%
if Cluster
    a = get_newaxes([xpos ypos(1,1) xwidth ypos(1,2)]);
    h = dendrogram(ClusterStats.Cltree,Nclust,'Reorder',...
        ClusterStats.ClleafOrder,'orientation','top');
    for i=1:length(h);set(h(i),'color','k'),end
    set(gca,'xtick',[],'ytick',[],'Visible','off')
    xlim([.5 Nbars+.5])

    txtoutput(end+(1:2),:) = [{'Dist with left Clust'; 'Dist with right Clust'} ...
        num2cellstr([0 diag(ClusterStats.ClDist,-1)' NaN;
        diag(ClusterStats.ClDist,1)' NaN NaN])];
end

GroupPos = 1:Nbars; % default

for iF=1:length(fields)
    a(iF+Cluster) = get_newaxes([xpos ypos(iF+Cluster,1) xwidth ypos(iF+Cluster,2)],1);

        % text output
        txtoutput{end+2,1} = fields{iF};

    if Fiscat(iF)
        h = bar(GroupPos, Ffraction{iF}, barw, 'stack');
        hl = legend(h(end:-1:1), Fcats{iF}(end:-1:1), ...
            'location','eastoutside','interpreter','none', plotopt.labels{:});
        posl = get(hl,'position');
        set(hl,'position',[xlpos posl(2:4)])

        if isfield(plotopt, [fields{iF} 'Colors']) && isfield(plotopt, fields{iF})
            for i=1:length(h)
                if any(strcmp(get(h(i),'DisplayName'), plotopt.(fields{iF})))
                    set(h(i),'facecolor',plotopt.([fields{iF} 'Colors'])( ...
                        strcmp(get(h(i),'DisplayName'), plotopt.(fields{iF})),:));
                else
                    set(h(i),'facecolor','w')
                end
            end
        end
        ylim([0 1.01])
        ylabel(['Dist. ' fields{iF}], 'interpreter','none', plotopt.labels{:})

        % text output
        txtoutput(end+(1:length(Fcats{iF})),:) = [Fcats{iF} num2cellstr(Ffraction{iF}')];
    else
        txtoutput{end,1} = [txtoutput{end,1} ' (quantiles)'];

        plot([0 max(Gidx)+1], [1 1], '-r')
        plot([0 max(Gidx)+1], [0 0], '-r')
        plot([0 max(Gidx)+1], -log10([.05 .05]), '-r')
        plot([0 max(Gidx)+1], [3 3], '-r')

        quants = [0 .1 .25 .5 .75 .9 1];
        idx = size(txtoutput,1)+(1:length(quants));
        txtoutput(idx,1) = num2cellstr(quants');

        for i=1:max(Gidx)
            plot_vbox(GroupPos(i), t_group_data.(fields{iF})(Gidx==i));
            txtoutput(idx,i+1) = ...
                num2cellstr(quantile(t_group_data.(fields{iF})(Gidx==i),quants));
        end
        ylims = [min(t_group_data.(fields{iF})) max(t_group_data.(fields{iF}))];
        ylim([ylims(1)-.05*diff(ylims) ylims(2)+.05*diff(ylims)])

        ylabel(fields{iF}, 'interpreter','none', plotopt.labels{:})
    end
    set(gca, plotopt.axis{:})

end

Clustercount = hist(Gidx,1:Nbars);
a(end+1) = get_newaxes([xpos ypos(end,1) xwidth ypos(end,2)],1);
bar(GroupPos, Clustercount, barw,'k');
ylim([0 max(Clustercount(1:Nclust))*1.05])


for i=1:length(a)
    xlim(a(i), [.5 Nbars+.5])
    set(a(i),'xtick',[])
end
set(a(end),'xtick',GroupPos,'xticklabel',1:Nbars, plotopt.axis{:})


%%

cell2tsv([plotopt.figurename '.tsv'], txtoutput);
