function plot_DrugDesing(Design, figurename, folder)
% plot_DrugDesing(Design, figurename, folder)

if ~exist('figurename','var')
    figurename = '';
end
if ~exist('folder','var')
    folder = '.';
end

for iD = 1:length(Design)

    get_newfigure(1e4+iD, [iD*50 20 700 600], [folder filesep ...
        sprintf('Design_%s_%i.pdf', figurename, iD)], ...
        'Name', sprintf('Design #%i', iD))

    nDisp = length(Design(iD).Drugs);
    if isfield(Design(iD), 'Perturbations')
        nDisp = nDisp + length(Design(iD).Perturbations);
    end

    nCols = ceil(sqrt(nDisp+1));
    nRows = ceil((nDisp+1)/nCols);

    xspacing = .03;
    axis_width = (1-(nCols+1)*xspacing)/nCols;
    yspacing = .1;
    axis_height = (1-(nRows+1)*yspacing)/nRows;

    colormap([.8 .1 .1; (1:-.01:0)'*[1 1 1]])

    for iDr = 1:length(Design(iD).Drugs)

        iC = mod(iDr-1, nCols)+1;
        iR = ceil(iDr/nCols);

        get_newaxes(axespos(iC, iR),1,...
            'fontsize',6);

        title(sprintf('%s [%.2g - %.2g]', Design(iD).Drugs(iDr).DrugName, ...
            min(Design(iD).Drugs(iDr).layout(Design(iD).Drugs(iDr).layout(:)>0)), ...
            max(Design(iD).Drugs(iDr).layout(:))), ...
            'fontsize',8,'fontweight','bold')

        xlim([.5 Design(iD).plate_dims(2)+.5])
        ylim([.5 Design(iD).plate_dims(1)+.5])

        set(gca,'ytick',1:2:16, 'yticklabel', cellfun(@(x) {char(x)},num2cell(65:2:80)),...
            'xtick',1:3:30,'ydir','reverse')

        if all(Design(iD).Drugs(iDr).layout(:)==0)
            continue
        end
        conc = log10(Design(iD).Drugs(iDr).layout);
        [range, conc] = DataRangeCap(conc);
        conc(~Design(iD).treated_wells) = NaN;

        imagesc(conc, range)

    end


    if isfield(Design(iD), 'Perturbations')


        for iPr = 1:length(Design(iD).Perturbations)

            iC = mod(iPr-1+length(Design(iD).Drugs), nCols)+1;
            iR = ceil((iPr+length(Design(iD).Drugs))/nCols);

            get_newaxes(axespos(iC, iR),1,...
                'fontsize',6);

            conc = Design(iD).Perturbations(iPr).layout;

            if islogical(conc)
                conc = conc+1;
                label = 'F/T';
            elseif iscell(conc) && iscellstr(conc(~cellfun(@isempty,conc)));
                conc(cellfun(@isempty,conc)) = {''};
                cats = unique(conc(:));
                val = NaN(size(conc));
                for i=1:length(cats)
                    val(strcmp(cats{i},conc)) = i-1;
                end
                conc = val/(length(cats)-.5);
                label = ['(' strjoin(cats,',') ')'];
            elseif (max(conc(:))/min(conc(:)))>10 && max(conc(:))*min(conc(:))>=0
                conc = log10(conc);
                label = sprintf('log10  [%.2g - %.2g]', max(conc(:)), min(conc(:)));
            else
                label = sprintf('[%.2g - %.2g]', max(conc(:)), min(conc(:)));
            end
            [range, conc] = DataRangeCap(conc);


            imagesc(conc, range)

            title({Design(iD).Perturbations(iPr).Name; label}, ...
                'fontsize',8,'fontweight','bold')

            xlim([.5 Design(iD).plate_dims(2)+.5])
            ylim([.5 Design(iD).plate_dims(1)+.5])

            set(gca,'ytick',1:2:16, 'yticklabel', cellfun(@(x) {char(x)},num2cell(65:2:80)),...
                'xtick',1:3:30,'ydir','reverse')
        end

    end

    get_newaxes(axespos(nCols, nRows),nCols*nRows,'fontsize',6);

    if ~isempty(Design(iD).Drugs)
        allDrugs = reshape([Design(iD).Drugs.layout], ...
            [Design(iD).plate_dims length(Design(iD).Drugs)]);
        ctrls = all(allDrugs==0,3)*1;
        ctrls(~Design(iD).treated_wells) = NaN;

        imagesc(ctrls, DataRangeCap([.2 1]))

        title('Control wells','fontsize',8,'fontweight','bold')

        xlim([.5 Design(iD).plate_dims(2)+.5])
        ylim([.5 Design(iD).plate_dims(1)+.5])

        set(gca,'ytick',1:2:16, 'yticklabel', cellfun(@(x) {char(x)},num2cell(65:2:80)),...
            'xtick',1:3:30,'ydir','reverse')

    end
end



    function pos = axespos(iC, iR)
        pos = [xspacing*1.5+(iC-1)*(axis_width+xspacing) ...
            yspacing*1+(nRows-iR)*(axis_height+yspacing) axis_width axis_height];
    end

    function [range, conc] = DataRangeCap(conc)
        range = [min(conc(~isinf(conc))) max(conc(~isinf(conc)))];
        if range(1)==range(2)
            range = range+[-.1 .1];
        end
        conc(conc==-Inf) = min(conc(conc>-Inf)) -.2*diff(range);
        range = [range(1)-.3*diff(range) range(2)+.1*diff(range)];
    end

end
