
function plot_DrugPairs_Bliss(t_mean)


t_DrugPairs = unique(t_mean(t_mean.DrugName2~='-',...
    {'CellLine' 'DrugName' 'DrugName2'}));

nRows = 2*floor(sqrt(height(t_DrugPairs)));
nCol = ceil(2*height(t_DrugPairs)/nRows);

global Plotting_parameters
for iDp=1:height(t_DrugPairs)

    t_pair = t_mean((eqtable(t_mean, t_DrugPairs(iDp, [2 3])) | ...
        (t_mean.DrugName==t_DrugPairs.DrugName(iDp) & t_mean.DrugName2=='-') | ...
        (t_mean.DrugName==t_DrugPairs.DrugName2(iDp) & t_mean.DrugName2=='-')) & ...
        t_mean.CellLine==t_DrugPairs.CellLine(iDp), :);

    [BlissScore, Bliss, Results, Concs] = EvaluateBliss(t_pair,[],1,.3);

    ridx = [1 find(any(~isnan(Bliss),2))'];
    cidx = [1 find(any(~isnan(Bliss)))];



    a = get_subaxes(nRows, nCol, 2*ceil(iDp/nCol)-1, mod1(iDp,nCol),0,...
        'yspacing',.05,'xspacing',.04,'xLshift',-.02);
    imagesc(Results(ridx,cidx), [-1.1 1.1])
    colormap(a, Plotting_parameters.cmapRB)

    title(char(t_DrugPairs.CellLine(iDp)), 'fontsize',8,'fontweight','bold', ...
        'interpreter','none')

    if iDp==1
        hc = colorbar;
        set(hc,'position',[.95 .62 .01 .3],Plotting_parameters.axes{:})
        ylabel(hc,'Relative growth', Plotting_parameters.axislabel{:})
        ylim(hc,[-1 1.1])
    end



    a(2) = get_subaxes(nRows, nCol, 2*ceil(iDp/nCol), mod1(iDp,nCol),0,...
        'yspacing',.05,'xspacing',.04,'xLshift',-.02);

    imagesc(Bliss(ridx,cidx), [-.4 .4])
    colormap(a(2), [.9 .9 .9; Plotting_parameters.cmapBBY])

    if iDp==1
        hc = colorbar;
        set(hc,'position',[.95 .12 .01 .3],Plotting_parameters.axes{:})
        ylabel(hc,'Bliss excess', Plotting_parameters.axislabel{:})
        ylim(hc,[-.2 .4])
    end

    xlabel(a(2), [char(t_DrugPairs.DrugName2(iDp)) ' (uM)'], Plotting_parameters.axislabel{:})
    for j=1:2
        ylabel(a(j), [char(t_DrugPairs.DrugName(iDp)) ' (uM)'], Plotting_parameters.axislabel{:})
        set(a(j), Plotting_parameters.axes{:}, 'xtick', 1:2:length(cidx), ...
            'xticklabel', [0; round(Concs{2}(cidx(3:2:end)-1),2)], ...
            'ytick', 1:2:length(ridx), ...
            'yticklabel', [0; round(Concs{1}(ridx(3:2:end)-1),2)])
    end
    set(a(2), 'position', get(a(2),'position')+[0 .02 0 0])
end
