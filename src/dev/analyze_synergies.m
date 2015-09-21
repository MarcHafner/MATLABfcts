function a = analyze_synergies(t_CL, t_Results, Drugs, removed_replicates, fignum, figtag)
%  a = analyze_synergies(t_CL, t_Results, Drugs, removed_replicates, figtag)
%
%
%   implemented only for IC cruves, bit GI


Generate_Plotting_parameters

if ~exist('removed_replicates','var') || isempty(removed_replicates)
    removed_replicates = [];
end

if ~exist('fignum','var') || isempty(fignum)
    fignum = 0;
end

if ~exist('figtag','var') || isempty(figtag)
    figtag = '';
end

assert(HasUniqueElement(t_CL.CellLine))
CellLine = char(unique(t_CL.CellLine));


%% plot the controls and the single drugs titration

removed_replicates = [];

CtrlIdx = t_CL.Ctrl==1;
meanCtrl = varfun(@mean,t_CL(CtrlIdx,:),'InputVariables', ...
    'Cellcount','GroupingVariables','Replicate');
meanCtrl = meanCtrl.mean_Cellcount;


if ismember('Day0', t_CL.Properties.VariableNames) && any(t_CL.Day0)
    DoGI50 = true;
    SeededNumber = trimmean(t_CL.Cellcount(t_CL.Day0==1),50);
else
    DoGI50 = false;
end

if isfield(Drugs,'Primary')
    PrimDrugCidx = find([Drugs.Primary]);
    SecDrugCidx = find(~[Drugs.Primary]);
    ComboDidx = Allpairs(PrimDrugCidx, SecDrugCidx);
    Nrows = length(PrimDrugCidx);
    Ncols = length(SecDrugCidx);
else
    ComboDidx = nchoosek(1:length(Drugs),2);
    Nrows = length(Drugs)-1;
    Ncols = length(Drugs)-1;
end

get_newfigure(fignum+1, [120 80 50+Ncols*150 40+170*Nrows], [CellLine '_RelCellCnt' figtag '.png'])
get_newfigure(fignum+3, [400 380 50+Ncols*150 40+170*Nrows], [CellLine '_Bliss' figtag '.png'])

for iD = 1:size(ComboDidx,1)
    Doses1 = unique([0 Drugs(ComboDidx(iD,1)).ComboDoses]);
    Drug1 = Drugs(ComboDidx(iD,1)).DrugName;
    Drug1CIdx = find(strcmp(t_CL.Properties.VariableNames, Drug1));

    Drug1Idx = table2array(t_CL(:,Drug1CIdx))~=0 ;

    Doses2 = unique([0 Drugs(ComboDidx(iD,2)).ComboDoses]);
    Drug2 = Drugs(ComboDidx(iD,2)).DrugName;

    Drug2CIdx = find(strcmp(t_CL.Properties.VariableNames, Drug2));


    ComboIdx = table2array(t_CL(:,Drug1CIdx))~=0 & ...
        table2array(t_CL(:,Drug2CIdx))~=0;

    %     Drug1fit = t_Results.ICfit{ComboDidx(iD,1)};
    %     Drug2fit = t_Results.ICfit{ComboDidx(iD,2)};

    Drug1fit = @(x) spline([0 t_Results.Doses{ComboDidx(iD,1)}], ...
        [1 smooth(t_Results.Rel_CellCnt{ComboDidx(iD,1)},'lowess')'], x);

    Drug2fit = @(x) spline([0 t_Results.Doses{ComboDidx(iD,2)}], ...
        [1 smooth(t_Results.Rel_CellCnt{ComboDidx(iD,2)},'lowess')'], x);

    CellCnt = NaN(length(Doses1), length(Doses2), max(t_CL.Replicate));

    AddRelcnt = NaN(length(Doses1), length(Doses2));

    for iDo1 = 1:length(Doses1)
        for iDo2 = 1:length(Doses2)
            if iDo1 == 1 && iDo2 == 1
                CellCnt(iDo1,iDo2,:) = meanCtrl;
                AddRelcnt(iDo1,iDo2) = 1;
            elseif iDo1 == 1
                CellCnt(iDo1,iDo2,:) = meanCtrl*Drug2fit(Doses2(iDo2));
                AddRelcnt(iDo1,iDo2) = min(1,Drug2fit(Doses2(iDo2)));
            elseif iDo2 == 1
                CellCnt(iDo1,iDo2,:) = meanCtrl*Drug1fit(Doses1(iDo1));
                AddRelcnt(iDo1,iDo2) = min(1,Drug1fit(Doses1(iDo1)));
            else
                AddRelcnt(iDo1,iDo2) = ...
                    min(1,Drug1fit(Doses1(iDo1)))*...
                    min(1,Drug2fit(Doses2(iDo2)));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%     in case of a flat fit, use the actual cell count
                %%%%     at the given point is lower
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for iR = setdiff(1:max(t_CL.Replicate),removed_replicates)
                    CellCnt(iDo1,iDo2,iR) = mean(t_CL.Cellcount(...
                        t_CL.(Drug1)==Doses1(iDo1) & ...
                        t_CL.(Drug2)==Doses2(iDo2) & ...
                        t_CL.Replicate==iR));
                end
            end
        end
    end

    Relcnt = nanmean(CellCnt./repmat(reshape(meanCtrl,1,1,[]),...
        length(Doses1),length(Doses2)),3);

    for iF=[1 3]
        figure(fignum+iF)
        get_newaxes([.08+.1/Ncols+.91*(ComboDidx(iD,2)-min(ComboDidx(:,2)))/Ncols ...
            .06+.12/Nrows+.85*(ComboDidx(iD,1)-min(ComboDidx(:,1)))/Nrows ...
            -.06+.85/Ncols -.09+.8/Nrows])

        if iF==1 % cell count
            clims = [.1 1.2];
            imagesc(Relcnt,clims)

            if DoGI50
                c0 = SeededNumber/mean(meanCtrl);
                l1 = clims(2)-c0;
                l2 = c0-clims(1);
                if l1>l2
                    idx = ceil(l2*length(Plotting_parameters.cmapWBr)/l1);
                    cmap = [Plotting_parameters.cmapWBr(idx:-1:1,:);
                        Plotting_parameters.cmapWP];
                else
                    idx = ceil(l1*length(Plotting_parameters.cmapWP)/l2);
                    cmap = [Plotting_parameters.cmapWBr(end:-1:1,:);
                        Plotting_parameters.cmapWP(1:idx,:)];
                end
                colormap(cmap)
            else
                colormap(Plotting_parameters.cmapWP)
            end


        elseif iF==2 % different predicted
            delta = -(min(Relcnt,1)-AddRelcnt);
            delta(abs(delta)<.05) = 0;

            imagesc(delta,[-.301 .301])

            colormap(Plotting_parameters.cmapBR)
        else    % Bliss combination index
            I1 = 1-AddRelcnt(1,:);
            I2 = 1-AddRelcnt(:,1);
            I12= 1-Relcnt;
            Bliss = I12-(ones(size(I2))*I1+I2*ones(size(I1))-I2*I1);
            Bliss(1,:) = 0;
            Bliss(:,1) = 0;

            imagesc(Bliss,[-.301 .301])

            colormap(Plotting_parameters.cmapBR)

        end
        if iD==1
            hc = colorbar('location','northoutside');
            set(hc,'position',[.4 .93-.1/Nrows .2 .02],'fontsize',8)
            if iF==1
                xlabel(hc,[CellLine ': Cell count (rel to ctrl)'],'fontsize',10,'fontweight','bold')
            elseif iF==2
                xlabel(hc,[CellLine ': Delta Rel cell count with prediction'],'fontsize',10,'fontweight','bold')
            elseif iF==3
                xlabel(hc,[CellLine ': Delta Bliss (red=synergy)'],'fontsize',10,'fontweight','bold')
            end
        end

        set(gca,'xtick',2:2:length(Doses2),'xticklabel',.1*round(10*log10(Doses2(2:2:end))), ...
            'ytick',2:length(Doses1),'yticklabel',.1*round(10*log10(Doses1(2:end))),'fontsize',8)
        ylabel([Drug1 ' log_{10}(\muM)'])
        xlabel([Drug2 ' log_{10}(\muM)'])
    end

end
