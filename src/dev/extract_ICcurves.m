% t_Results = extract_ICcurves(t_CL, Drugs, removed_replicates, fignum, figtag, opt)
%   opt.['pcutoff' 'plotGI50']
%
function [t_Results, a] = extract_ICcurves(t_CL, Drugs, removed_replicates, fignum, figtag, opt)


if ~exist('removed_replicates','var') || isempty(removed_replicates)
    removed_replicates = [];
end


if ~exist('fignum','var') || isempty(fignum)
    fignum = 0;
end

if ~exist('figtag','var') || isempty(figtag)
    figtag = '';
end


pcutoff = .05;

if exist('opt','var')
    fields = {'pcutoff' 'plotGI50'};
    for field = fields
        if isfield(opt,field{:})
            eval([field{:} ' = opt.' field{:} ';'])
        end
    end
end

fitopt.pcutoff = pcutoff;

%%
assert(HasUniqueElement(t_CL.CellLine))
CellLine = char(unique(t_CL.CellLine));

%% plot the controls and the single drugs titration


CtrlIdx = t_CL.Ctrl==1;

GoodReplicates = ToRow(unique(t_CL.Replicate));
if ~isempty(removed_replicates)
    GoodReplicates = setdiff(ToRow(unique(t_CL.Replicate)), removed_replicates);
end
GoodReplicates = GoodReplicates(GoodReplicates>0);

if ismember('Day0', t_CL.Properties.VariableNames) && any(t_CL.Day0)
    DoGI50 = true;
    Day0Cnt = trimmean(t_CL.Cellcount(t_CL.Day0==1),50);
else
    DoGI50 = false;
end
if ~exist('plotGI50','var')
    plotGI50 = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     remove the edges if there are made of controls and there are
%%%         enough controls in the center
%%%%%%%%%%%%%%%%

if fignum
    get_newfigure(fignum, [50 250 900 500],  ...
        [CellLine '_ICcurves' figtag '.pdf'])
end
fprintf('Fit for %s:\n',CellLine);
ctrls = NaN(1,max(t_CL.Replicate));
for iR = 1:max(t_CL.Replicate)
    RepIdx = t_CL.Replicate==iR;

    ctrls(iR) = trimmean(t_CL.Cellcount(CtrlIdx & RepIdx),50);
    std_ctrl = std(t_CL.Cellcount(CtrlIdx & RepIdx));

    fprintf('\tUntreated ctrl (rep %i): %.0f +/- %.0f (%i)\n', ...
        iR, ctrls(iR), std_ctrl,sum(CtrlIdx & RepIdx));

end



if DoGI50
    fprintf('\tSeeding number: %.0f +/- %.0f\n', Day0Cnt, ...
        std(t_CL.Cellcount(t_CL.Day0==1)));
    GI50 = NaN(1, length(Drugs));
    if fignum && plotGI50
        get_newfigure(fignum+1, [50 250 900 500],  ...
            [CellLine '_GIcurves_' figtag '.pdf'])
    end
    t_GI = table;
end

t_IC = table;
t_labels = table;

for iD=1:length(Drugs)
    DrugName = Drugs(iD).DrugName;
    t_labels = [t_labels;
        table({CellLine}, {DrugName}, 'VariableNames', {'CellLine' 'Drug'})];

    DrugCidx = strcmp(t_CL.Properties.VariableNames,DrugName);
    otherDrugCidx = ismember(t_CL.Properties.VariableNames,...
        {Drugs(setdiff(1:length(Drugs),iD)).DrugName});

    DrugIdx = table2array(t_CL(:,DrugCidx))~=0 & ...
        all(table2array(t_CL(:,otherDrugCidx))==0,2);

    Doses = unique(t_CL.(DrugName)(t_CL.Replicate==1 & DrugIdx));
    assert(all(Doses==unique(setdiff(t_CL.(DrugName)(DrugIdx),0))))

    Relcnt = NaN(max(t_CL.Replicate), length(Doses));

    if DoGI50
        Relgrowth = NaN(max(t_CL.Replicate), length(Doses));
    end
    if fignum
        xlims = [min(log10(Doses))-.1 max(log10(Doses))+.1]*1.05;

        figure(fignum)
        a(iD) = subplot(floor(sqrt(length(Drugs))), ...
            ceil(length(Drugs)/floor(sqrt(length(Drugs)))),iD);
        plot(xlims, [1 1],'-c')
        hold on

        if plotGI50
            figure(fignum+1)
            subplot(floor(sqrt(length(Drugs))), ...
                ceil(length(Drugs)/floor(sqrt(length(Drugs)))),iD)
            plot(xlims, [1 1],'-c')
            hold on
        end
    end

    for iR = 1:max(t_CL.Replicate)
        RepIdx = t_CL.Replicate==iR;

        ctrl = mean(t_CL.Cellcount(CtrlIdx & RepIdx));
        std_ctrl = std(t_CL.Cellcount(CtrlIdx & RepIdx));

        t_doses = sortrows(t_CL(RepIdx & DrugIdx,:),DrugName);
        assert(all(unique(t_doses.(DrugName))==Doses))

        Relcnt(iR,:) = varfun(@mean,t_doses,'GroupingVar',DrugName, ...
            'InputVar', 'Cellcount', 'outputformat','uniform')/ctrls(iR);
        if fignum
            figure(fignum)
            if ismember(iR, GoodReplicates)
                plot(log10(Doses),Relcnt(iR,:),'.-');
            else
                plot(log10(Doses),Relcnt(iR,:),'.-','color',[.7 .7 .7]);
            end

            plot(xlims, [1 1]-(std_ctrl/ctrl),':c')
            plot(xlims, [1 1]+(std_ctrl/ctrl),':c')
        end

        if DoGI50
            Relgrowth(iR,:) = (varfun(@mean,t_doses,'GroupingVar',DrugName,...
                'InputVar', 'Cellcount', 'outputformat','uniform')-Day0Cnt)/...
                (ctrls(iR)-Day0Cnt);

            if fignum && plotGI50
                figure(fignum+1)
                if ismember(iR, GoodReplicates)
                    plot(log10(Doses),Relgrowth(iR,:),'.-');
                else
                    plot(log10(Doses),Relgrowth(iR,:),'.-','color',[.7 .7 .7]);
                end

                plot(xlims, [1 1]-(std_ctrl/ctrl),':c')
                plot(xlims, [1 1]+(std_ctrl/ctrl),':c')
            end
        end
    end


    [IC50, Hill, ~, Emax, Area, r2, EC50, fit, ~, log] = ...
        ICcurve_fit(Doses, mean(Relcnt(GoodReplicates,:)), 'IC50', fitopt);

    t_IC = [t_IC;
        cell2table({IC50, Hill, Emax, Area, r2, EC50, fit, log, {Doses'}, ...
        {mean(Relcnt(GoodReplicates,:))}, mean(ctrls(GoodReplicates))}, 'VariableNames', ...
        {'IC50' 'Hill' 'Emax' 'Area' 'r2' 'EC50' 'ICfit' 'log' 'Doses' 'Rel_CellCnt' 'CtrlCnt'})];


    if DoGI50
        [GI50, Hill, ~, GI_Emax, Area, r2, ~, GI_fit, ~, log] = ...
            ICcurve_fit(Doses, mean(Relgrowth(GoodReplicates,:)), 'GI50', fitopt);

        t_GI = [t_GI;
            cell2table({GI50, Hill, GI_Emax, Area, r2, GI_fit, log, ...
            {mean(Relgrowth(GoodReplicates,:))} Day0Cnt*ones(size(r2))},...
            'VariableNames', ...
            {'GI50' 'Hill_GI' 'Emax_GI' 'Area_GI' 'r2_GI' 'GIfit' ...
            'log_GI' 'Rel_growth' 'Day0Cnt'})];
    end

    if fignum
        figure(fignum)
        plot(log10(Doses), mean(Relcnt(GoodReplicates,:)), '.-k')
        plot(log10(Doses), fit(Doses), '.-r','linewidth',2)
        plot(log10(IC50)*[1 1], [0 .5*(1+Emax)], '-r')
        plot(xlims, [1 1]*Emax, '-r')
        score = sprintf('log_{10}(IC_{50})=%.2g', log10(IC50));

        if DoGI50
            plot(xlims, [1 1]*Day0Cnt/ctrl,':k')
            plot(log10(GI50)*[1 1], [0 .5*(1+Day0Cnt/ctrl)], '-k')
            score = [score sprintf(', log_{10}(GI_{50})=%.2g ', log10(GI50))];
        end

        title({CellLine; DrugName; score},...
            'fontsize',8,'fontweight','bold')
        set(gca,'ytick',0:.2:1, 'fontsize',6,'xtick',...
            (.1*floor(10*min(log10(Doses))):.5:(.1*ceil(10*max(log10(Doses))))))
        xlabel('log_{10}(Dose) (\muM)','fontsize',6,'fontweight','bold')
        ylabel('Relative cell cnt','fontsize',6,'fontweight','bold')
        ylim([0 1.4])
        xlim(xlims)
    end

    if fignum && plotGI50
        figure(fignum+1)
        plot(log10(Doses), mean(Relgrowth(GoodReplicates,:)), '.-k')
        plot(log10(Doses), GI_fit(Doses), '.-r','linewidth',2)


        title({CellLine; DrugName; sprintf('r=%.2f', r2)},...
            'fontsize',8,'fontweight','bold')
        set(gca,'ytick',-2:.25:1, 'fontsize',6,'xtick',...
            (.1*floor(10*min(log10(Doses))):.5:(.1*ceil(10*max(log10(Doses))))))
        xlabel('log_{10}(Dose) (\muM)','fontsize',6,'fontweight','bold')
        ylabel('Relative growth','fontsize',6,'fontweight','bold')
        ylim([min([-.1 min(Relgrowth(GoodReplicates,:))]) 1.4])
        xlim(xlims)

    end
end


if DoGI50
    t_Results = [t_labels t_IC t_GI];
else
    t_Results = [t_labels t_IC];
end
