function t_CLcorrected = correct_edges(t_CL, figtag)
% t_CLcorrected = correct_edges(t_CL, figtag)
%
%


if ~exist('figtag','var')
    figtag = 0;
end

assert(HasUniqueElement(t_CL.CellLine))
CellLine = char(unique(t_CL.CellLine));

t_CLcorrected = t_CL;

fprintf('Correction of edges for %s\n', CellLine)

bins = [min(t_CL.Cellcount) max(t_CL.Cellcount)];
bins = [floor(.95*bins(1)):floor(diff(bins)/30):ceil(1.05*bins(2))];


if figtag
    get_newfigure(figtag,[100 100 700 450])
end
ratios = NaN(1,max(t_CL.Replicate));
pvals = ratios;
nR = max(t_CL.Replicate);
for iR = 1:nR
    t_ctrl_plate = t_CL(t_CL.Replicate==iR & t_CL.Ctrl==1,...
        {'Row','Column','Cellcount'});

    if figtag
        a1 = get_newaxes([.04+(iR-1)/nR .4 -.05+.95/nR .5],1);
        a2 = get_newaxes([.04+(iR-1)/nR .05 -.05+.95/nR .3],1);
    else
        a1 =[];a2 = [];
    end

    [ratios(iR), pvals(iR)] = compare_edge_center(t_ctrl_plate, bins, a1, a2, CellLine);

    idxes = t_CL.Replicate==iR;
    t_CLcorrected(idxes,:) = apply_edge_correction(t_CL(idxes,:), ...
        ratios(iR), pvals(iR));
end


if figtag
    get_newfigure(1+figtag,[100 500 700 450])
end
nD0 = max(t_CL.Day0);
for iD0 = 1:max(t_CL.Day0)
    idxes = t_CL.Day0==iD0;
    t_ctrl_plate = t_CL(idxes,{'Row','Column','Cellcount'});

    if figtag
        a1 = get_newaxes([.04+(iD0-1)/nD0 .4 -.05+.95/nD0 .5],1);
        a2 = get_newaxes([.04+(iD0-1)/nD0 .05 -.05+.95/nD0 .3],1);
    else
        a1 =[];a2 = [];
    end

    [ratio, pval] = compare_edge_center(t_ctrl_plate, bins, a1, a2, CellLine);

    t_CLcorrected(idxes,:) = apply_edge_correction(t_CL(idxes,:), ...
        ratio, pval);
end

if figtag
    get_newfigure(2+figtag,[100 500 700 450])
end
nUnt = max(t_CL.Untrt);
for iUnt = 1:nUnt
    idxes = t_CL.Untrt==iUnt;
    t_ctrl_plate = t_CL(idxes,{'Row','Column','Cellcount'});

    if figtag
        a1 = get_newaxes([.04+(iUnt-1)/nUnt .4 -.05+.95/nUnt .5],1);
        a2 = get_newaxes([.04+(iUnt-1)/nUnt .05 -.05+.95/nUnt .3],1);
    else
        a1 =[];a2 = [];
    end

    [ratio, pval] = compare_edge_center(t_ctrl_plate, bins, a1, a2, CellLine);
    %%% TO IMPROVE: IMPLEMENT A CONTROL FOR ROW BIAS --> FLAG (NO
    %%% CORRECTION

    t_CLcorrected(idxes,:) = apply_edge_correction(t_CL(idxes,:), ...
        ratio, pval);
end


end

function [ratio, pval] = compare_edge_center(t_ctrl_plate, bins, a1, a2, CellLine)

Generate_Plotting_parameters

edgeCtrl = t_ctrl_plate( t_ctrl_plate.Row==1 | t_ctrl_plate.Row==16 | ...
    t_ctrl_plate.Column==1 | t_ctrl_plate.Column==24, :);
if isempty(edgeCtrl)
    return
end

centerCtrl = setdiff(t_ctrl_plate,edgeCtrl);

pval = ranksum(centerCtrl.Cellcount,edgeCtrl.Cellcount);
ratio =  median(centerCtrl.Cellcount)/median(edgeCtrl.Cellcount);

if isempty(a1) || isempty(a2)
    return
end

set(gcf,'CurrentAxes',a1)

n = ksdensity(edgeCtrl.Cellcount,bins,'width',diff(bins([1 2]))*2);
h = plot(bins,n);

n = ksdensity(centerCtrl.Cellcount,bins,'width',diff(bins([1 2]))*2);
h(2) = plot(bins,n,'r');

legend(h, {['Edge, n=' num2str(height(edgeCtrl))] ...
    ['Center, n=' num2str(height(centerCtrl))]})
xlim([min(bins) max(bins)])
set(gca,'ytick',[])
title({ sprintf('%s; CV=%.2f, center cnt=%.0f',CellLine, ...
    std(t_ctrl_plate.Cellcount)/mean(t_ctrl_plate.Cellcount), median(centerCtrl.Cellcount));
    sprintf('\\Delta=%.1f, ratio=%.3f, p=%.3f', ...
    median(centerCtrl.Cellcount)-median(edgeCtrl.Cellcount), ...
    ratio, pval)});

missingCol = setdiff(1:24, t_ctrl_plate.Column);
missingRow = setdiff(1:16, t_ctrl_plate.Row);
t_ctrl_plate = [t_ctrl_plate; cell2table(num2cell(...
    [ones(length(missingCol),1) missingCol' NaN(length(missingCol),1);
    missingRow' ones(length(missingRow),1) NaN(length(missingRow),1)]), ...
    'VariableNames', t_ctrl_plate.Properties.VariableNames)];


set(gca,'fontsize',8)


[Ctrl2D, labels] = table_to_ndarray(t_ctrl_plate,...
    'KeyVars', {'Row' 'Column'}, 'ValVars', {'Cellcount'}, 'outer',1);
order = sortidx(labels{2}.Column);

set(gcf,'CurrentAxes',a2,'color','w')

imagesc(Ctrl2D(:,order), [min(bins) max(bins)])
xlim([.4 height(labels{2})+.6])
ylim([.4 height(labels{1})+.6])

colormap(Plotting_parameters.cmapWP)

set(gca,'fontsize',8,'xtick',[1 5 10 15 20 24],'ytick',[1 4 8 12 16],...
    'yticklabel',{'A' 'D' 'H' 'L' 'P'},'ydir','reverse','box','on');

end


function  subt_CLcorrected = apply_edge_correction(subt_CL, ratio, pval)

subt_CLcorrected = subt_CL;

if pval<.1
    edgeIdx = (subt_CL.Row==1 | subt_CL.Row==16 | ...
        subt_CL.Column==1 | subt_CL.Column==24) ;
    subt_CLcorrected.Cellcount(edgeIdx) = ...
        subt_CL.Cellcount(edgeIdx)*ratio;

    fprintf('\tplate corrected (p=%.2f, ratio=%.2f)\n', ...
        pval, ratio);
else
    fprintf('\tplate has NOT been corrected (p=%.2f, ratio=%.2f)\n', ...
        pval, ratio);
end
end
