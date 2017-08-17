function [LiveCells, DeadCells, Gates, AliveIdx, LDRlims, DNAlims] = DeadCellFilter(LDRtxt, varargin)

p = inputParser;

addOptional(p, 'DNA', [], @(x) isvector(x) & length(x)==length(LDRtxt))
addParameter(p, 'plotting', false, @islogical)
addParameter(p, 'interactive', false, @islogical)
addParameter(p, 'xDNA', 2.5:.02:8, @isvector)
addParameter(p, 'xLDR', -.01:.0002:(max(LDRtxt)+.01), @isvector)
addParameter(p, 'LDRcutoff', [], @isscalar)
addParameter(p, 'nsmooth', 5, @isnumeric)
addParameter(p, 'DNApks', [NaN NaN], @(x) isvector(x) & length(x)==2)
addParameter(p, 'Gates', NaN(2), @(x) ismatrix(x) & all(size(x)==2) & all(~isnan(x(:,1))))
addParameter(p, 'savefigure', '', @ischar)
addParameter(p, 'LDRlims', [], @(x) all(size(x)==[1 2]) && x(2)>x(1))
addParameter(p, 'DNAlims', [], @(x) all(size(x)==[1 2]) && x(2)>x(1))

parse(p,varargin{:});
p = p.Results;

if p.interactive, p.plotting = true; end
if ~isempty(p.savefigure), p.plotting = true; end

DNA = p.DNA;
p = rmfield(p,'DNA');
useDNA = ~isempty(DNA);

%%

% start with the LDR channel
f = ksdensity(LDRtxt, p.xLDR, 'width', 2.5*diff(p.xLDR(1:2)));

if ~isnan(p.Gates(1,2))
    % already defined gate
    Gates = p.Gates(1,:);
    if isnan(p.Gates(1,1)), Gates(1,1) = -Inf; end
else
    % determine the spread of the LDRtxt data to define cutoff
    [~, pk, LDRwdth] = findpeaks(f,'npeaks',1,'widthreference','halfprom','sortstr','descend');
    
    [~, minpk] = findpeaks(-f(pk:end),'npeaks',1); minpk=minpk+pk-1;
    LDRcutoff = p.xLDR(ceil(max(min(minpk, pk+6*LDRwdth), pk+2*LDRwdth)));
    Gates = [-Inf LDRcutoff];
end


if p.plotting
    % plotting results
    currfig = gcf;
    
    plot_pos = [
    .07 .6 .4 .37;
    .57 .6 .4 .37;
    .15 .05 .4 .4
    .6 .15 .3 .3];

    get_newfigure(45674,[505 600 550 500])
    
    if ~isempty(p.LDRlims), LDRlims = p.LDRlims; else
    LDRlims = quantile(LDRtxt, [5e-3 .995])+[-1 1]*2.5*diff(p.xLDR(1:2)); end
    
    % plot original data
    get_newaxes(plot_pos(1,:),1)
    plot(p.xLDR, log10(f+max(f)/100)-log10(max(f)/100))
    
    plot([Gates(1,2) max(Gates(1,1), min(p.xLDR))*[1 1] [1 1]*Gates(1,2)], ...
        [0 0 .5 .5 0]*log10(max(f)), '-', 'color', [.6 .9 .1]);
    pltgt1 = plot([Gates(1,2) max(Gates(1,1), min(p.xLDR))*[1 1] [1 1]*Gates(1,2)], ...
        [0 0 .5 .5 0]*log10(max(f)), 'r-');
    
    xlim(LDRlims)
    ylim([0 log10(max(f))-log10(max(f)/100)+.1])
    set(gca,'xtick',[],'ytick',[])
    
    pieax = get_newaxes(plot_pos(4,:));   
    set(gca,'xtick',[],'ytick',[],'visible','off') 
end


if useDNA
    % work with the DNA content if provided
    
    % log10 domain and capping
    logDNA = log10(min(max(DNA, 10^p.xDNA(3)),10^p.xDNA(end-2)));
    f2 = ksdensity(logDNA,p.xDNA);
    
    if p.plotting
        % plot original data
        get_newaxes(plot_pos(2,:),1)
        if ~isempty(p.DNAlims), DNAlims = p.DNAlims; else
        DNAlims = quantile(logDNA, [5e-3 .995])+[-1 1]*2.5*diff(p.xDNA(1:2)); end
        plot(p.xDNA, f2, '-k')
        xlim(DNAlims)
        set(gca,'xtick',[],'ytick',[])
    end
    
    
    % take only the cells with low LDR
    f2s = ksdensity(logDNA(LDRtxt>=Gates(1,1) & LDRtxt<=Gates(1,2)),p.xDNA);
    if p.plotting, plot(p.xDNA, f2s, '--r'), end
    
    [pks, idx] = findpeaks(f2s, 'sortstr', 'descend');
    idx = idx(pks>max(pks/10)); % remove lesser peaks
    DNAPks = p.xDNA(idx(1:min(3,length(idx)))); % take the 3 highest peaks with low EdU
    
    % find the DNA peak for G1
    if length(DNAPks)>1
        % more than one candidate
        if ~isnan(p.DNApks(1))
            % input matching G1 peak
            DNAPks = DNAPks(argmin(abs(DNAPks-p.DNApks(1))));
        elseif ~isnan(p.DNApks(2))
            % input matching S peak
            DNAPks = max(DNAPks(DNAPks<p.DNApks(2)));
        else
            % take the highest peak (most likely case in doubt)
            DNAPks = DNAPks(1);
        end
    end
    
    
    
    % get the cells in G2 phase
    hD = logDNA>DNAPks(1)+.4*log10(2) & LDRtxt>=Gates(1,1) & LDRtxt<=Gates(1,2);
    if any(hD)
        % found some cells in G2 phase
        f3 = ksdensity(logDNA(hD),p.xDNA);
        
        if p.plotting
            plot(p.xDNA, f3, ':')
        end
        
        [pks, idx] = findpeaks(smooth(f3,p.nsmooth), 'sortstr', 'descend');
        idx = idx(pks>max(pks/10)); % remove lesser peaks
        hDNAPks = p.xDNA(idx);
        % should be around log10(2) above the DNA peak in G1
        hDNAPks = hDNAPks(hDNAPks>(DNAPks(1)+.5*log10(2)));
        
        if length(hDNAPks)>1
            % more than one candidate
            if ~isnan(p.DNApks(2))
                % 2D analysis matching G2 peak
                hDNAPks = hDNAPks(argmin(abs(hDNAPks-PhasesCandidates(3,1))));
            else
                % take the peak closest to a 2-fold(most likely case in doubt)
                hDNAPks = hDNAPks(argmin(abs(hDNAPks-DNAPks(1)-log10(2))));
            end
        end
        
        DNAPks = [DNAPks hDNAPks];
        
    else
        % no cells found in G2 phase
        % default value (2-fold)
        DNAPks = [DNAPks DNAPks(1)+log10(2)];
        
    end
    
    
    if any(isnan(p.Gates(2,:)))
        % define areas
        Gates(2,:) = DNAPks + [-1.2 1.2]*diff(DNAPks);
    else
        Gates(2,:) = p.Gates(2,:);
    end
    
    
    if p.plotting
        plot(Gates(2,[1 1 2 2]), [0 max(f2)*[1.02 1.02] 0], '--', 'color', [.6 .9 .1]);
        pltgt2 = plot(Gates(2,[1 1 2 2]), [0 max(f2)*[1.02 1.02] 0], '-r');
        phases = {'G1'  'G2'};
        for i=1:2
            text(DNAPks(i), max(f2)*1.1, phases{i}, ...
                'fontsize', 14, 'fontweight', 'bold', 'color', [.6 .9 .1], ...
                'horizontalalign','center')
        end
        ylim([0 max(f2)*1.2])
        % plots with both channels
        
        get_newaxes(plot_pos(3,:),1)
        dscatter(logDNA, LDRtxt, 'MSIZE', 15, 'marker', 'o')
        
        plot(Gates(2,[1 1 2 2 1]), max(Gates(1,[1 2 2 1 1]),0), '--', ...
            'color', [.6 .9 .1], 'linewidth', 2);
        pltgt3 = plot(Gates(2,[1 1 2 2 1]), max(Gates(1,[1 2 2 1 1]),0), ...
            '-r', 'linewidth', 2);
        plot(DNAPks, [0 0], 'xk')
        plot(DNAPks, [0 0], 'ok', 'markersize', 14)
        xlim(DNAlims)
        ylim(LDRlims)
    end
    
    
end

% finalize the results
[LiveCells, DeadCells, AliveIdx] = EvalAliveIdx();



if p.interactive
    figpos = get(gcf,'position');
    
    minLDR = uicontrol('style', 'slider', 'callback', {@setGates,1});
    minLDR.Units = 'normalized';
    minLDR.InnerPosition = [plot_pos(1,1)-15/figpos(3) plot_pos(1,2)-.04 plot_pos(1,3)/2+30/figpos(3) .03];
    minLDR.Value = max(0,(Gates(1,1)-LDRlims(1))/diff(LDRlims));
    
    maxLDR = uicontrol('style', 'slider', 'callback', {@setGates,2});
    maxLDR.Units = 'normalized';
    maxLDR.Position = [plot_pos(1,1)-15/figpos(3) plot_pos(1,2)-.08 plot_pos(1,3)+30/figpos(3) .03];
    maxLDR.Value = (Gates(1,2)-LDRlims(1))/diff(LDRlims);
    
    if useDNA
        minDNA = uicontrol('style', 'slider', 'callback', {@setGates,3});
        minDNA.Units = 'normalized';
        minDNA.Position = [plot_pos(2,1)-15/figpos(3) plot_pos(2,2)-.04 plot_pos(2,3)+30/figpos(3) .03];
        minDNA.Value = (Gates(2,1)-DNAlims(1))/diff(DNAlims);
        
        maxDNA = uicontrol('style', 'slider', 'callback', {@setGates,4});
        maxDNA.Units = 'normalized';
        maxDNA.Position = [plot_pos(2,1)-15/figpos(3) plot_pos(2,2)-.08 plot_pos(2,3)+30/figpos(3) .03];
        maxDNA.Value = (Gates(2,2)-DNAlims(1))/diff(DNAlims);
    end
    
    approve = uicontrol('style', 'pushbutton');
    approve.Units = 'normalized';
    approve.Position = [.65 .03 .3 .1];
    approve.String = 'Approve';
    approve.Callback = @approveGate;
    
    waitfor(approve, 'backgroundcolor', 'g')
    
end

if p.plotting
    if ~isempty(p.savefigure)
        set(gcf,'Renderer','painters')
        saveas(gcf,p.savefigure)
    end
    
    figure(currfig)
end
%%
    function setGates(src, event, x)
        if x<3
            Gates(1,x) = (diff(LDRlims)*src.Value)+LDRlims(1);
        else
            Gates(2,x-2) = (diff(DNAlims)*src.Value)+DNAlims(1);
        end
        set(pltgt1, 'XData', [Gates(1,2) max(Gates(1,1), min(p.xLDR))*[1 1] [1 1]*Gates(1,2)])
        if useDNA
            set(pltgt2, 'XData', Gates(2,[1 1 2 2]))
            set(pltgt3, 'XData', Gates(2,[1 1 2 2 1]), 'YData', max(Gates(1,[1 2 2 1 1]),0));
        end
        
        [LiveCells, DeadCells, AliveIdx] = EvalAliveIdx();
    end

    function [alive, dead, idx] = EvalAliveIdx()
        idx = LDRtxt>=Gates(1,1) & LDRtxt<=Gates(1,2);
        if useDNA
            idx = idx & (logDNA>=Gates(2,1) & logDNA<=Gates(2,2));
        end
        alive = sum(idx);
        dead = sum(~idx);
        
        if p.plotting
            set(gcf,'currentaxes', pieax)
            cla
            ptxt = pie([alive dead] ,{'Live cells' sprintf('Dead cells (%.0f%%)', ...
                100*dead/(alive+dead))});
            set(ptxt(4),'fontsize',12, 'fontweight','bold')
        end
    end

    function approveGate(src, event)
        minLDR.Visible = 'off';
        maxLDR.Visible = 'off';
        if useDNA
            minDNA.Visible = 'off';
            maxDNA.Visible = 'off';
        end
        set(src, 'backgroundcolor', 'g')
    end


end