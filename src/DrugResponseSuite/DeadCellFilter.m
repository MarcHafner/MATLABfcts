function [LiveCells, DeadCells, Gates, AliveIdx] = DeadCellFilter(LDRtxt, varargin)
currfig = gcf;

p = inputParser;
addOptional(p, 'DNA', [], @(x) isvector(x) & length(x)==length(LDRtxt))
addParameter(p, 'plotting', false, @islogical)
addParameter(p, 'xDNA', 2.5:.02:8, @isvector)
addParameter(p, 'xLDR', -.01:.0002:(max(LDRtxt)+.01), @isvector)
addParameter(p, 'LDRcutoff', [], @isscalar)
addParameter(p, 'nsmooth', 5, @isnumeric)
addParameter(p, 'DNApks', [NaN NaN], @(x) isvector(x) & length(x)==2)
addParameter(p, 'Gates', NaN(2), @(x) ismatrix(x) & all(size(x)==2) & all(~isnan(x(:,1))))

parse(p,varargin{:});
p = p.Results;
DNA = p.DNA;
p = rmfield(p,'DNA');
useDNA = ~isempty(DNA);
%%

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

plot_pos = [
    .05 .6 .35 .37;
    .6 .6 .35 .37;
    .05 .1 .35 .37];

%
if p.plotting
    get_newfigure(45674,[505 600 500 650])
    ylims = quantile(LDRtxt, [1e-3 .999]);
    
    % plot original data
    get_newaxes(plot_pos(1,:),1)
    plot(p.xLDR, log10(f+max(f)/100)-log10(max(f)/100))
    pltgt1 = plot([Gates(1,2) max(Gates(1,1), min(p.xLDR))*[1 1] [1 1]*Gates(1,2)], ...
        [0 0 .5 .5 0]*log10(max(f)), 'r-');
    
    xlim(ylims)
    ylim([0 log10(max(f))-log10(max(f)/100)+.1])
    
end

AliveIdx = LDRtxt>=Gates(1,1) & LDRtxt<=Gates(1,2);

if useDNA
    % work with the DNA content
    
    % log10 domain and capping
    logDNA = log10(min(max(DNA, 10^p.xDNA(3)),10^p.xDNA(end-2)));
    f2 = ksdensity(logDNA,p.xDNA);
    
    if p.plotting
        % plot original data
        get_newaxes(plot_pos(2,:),1)
        hold on
        xlims = quantile(logDNA, [1e-3 .999])+[-.1 .2];
        plot(p.xDNA, f2, '-k')
        xlim(xlims)
        
    end
    
    
    % take only the cells with low LDR
    f2s = ksdensity(logDNA(AliveIdx),p.xDNA);
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
    
    
    
    % get the G2
    hD = logDNA>DNAPks(1)+.4*log10(2) & AliveIdx;
    if any(hD)
        % found some G2 cells
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
        % no G2 cells
        DNAPks = [DNAPks DNAPks(1)+log10(2)];
        
    end
    
    
    if any(isnan(p.Gates(2,:)))
        % define areas
        Gates(2,:) = DNAPks + [-1.2 1.2]*diff(DNAPks);
        AliveIdx = AliveIdx & (logDNA>=Gates(2,1) & logDNA<=Gates(2,2));
    else
        Gates(2,:) = p.Gates(2,:);
        AliveIdx = AliveIdx & (logDNA>=Gates(2,1) & logDNA<=Gates(2,2));
    end
    
    
    if p.plotting
        pltgt2 = plot(Gates(2,[1 1 2 2]), [0 max(f2)*[1.02 1.02] 0], '-r');
        phases = {'G1'  'G2'};
        for i=1:2
            text(DNAPks(i), max(f2)*1.1, phases{i}, ...
                'fontsize', 14, 'fontweight', 'bold', 'color', 'r', ...
                'horizontalalign','center')
        end
        ylim([0 max(f2)*1.2])
        % plots with both channels
        
        get_newaxes(plot_pos(3,:),1)
        dscatter(logDNA, LDRtxt, 'MSIZE', 15, 'marker', 'o')
        hold on
        ylim(ylims)
        xlim(xlims)
        
        pltgt3 = plot(Gates(2,[1 1 2 2 1]), max(Gates(1,[1 2 2 1 1]),0), '-r');
        plot(DNAPks, [0 0], 'xk')
        plot(DNAPks, [0 0], 'ok', 'markersize', 14)
        xlim(xlims)
        ylim(ylims)
    end
    
    
end
% finalize the results

LiveCells = sum(AliveIdx);
DeadCells = sum(~AliveIdx);


if p.plotting
    figpos = get(gcf,'position');
    
    minLDR = uicontrol('style', 'slider', 'callback', {@setGates,1});
    minLDR.Units = 'normalized';
    minLDR.InnerPosition = [plot_pos(1,1)-15/figpos(3) plot_pos(1,2)-.03 plot_pos(1,3)/2+30/figpos(3) .02];
    minLDR.Value = max(0,(Gates(1,1)-ylims(1))/diff(ylims));
    minLDR.Visible = 'off';
    
    maxLDR = uicontrol('style', 'slider', 'callback', {@setGates,2});
    maxLDR.Units = 'normalized';
    maxLDR.Position = [plot_pos(1,1)-15/figpos(3) plot_pos(1,2)-.06 plot_pos(1,3)+30/figpos(3) .02];
    maxLDR.Value = (Gates(1,2)-ylims(1))/diff(ylims);
    maxLDR.Visible = 'off';
    
    if useDNA
        minDNA = uicontrol('style', 'slider', 'callback', {@setGates,3});
        minDNA.Units = 'normalized';
        minDNA.Position = [plot_pos(2,1)-15/figpos(3) plot_pos(2,2)-.03 plot_pos(2,3)+30/figpos(3) .02];
        minDNA.Value = (Gates(2,1)-xlims(1))/diff(xlims);
        minDNA.Visible = 'off';
        
        maxDNA = uicontrol('style', 'slider', 'callback', {@setGates,4});
        maxDNA.Units = 'normalized';
        maxDNA.Position = [plot_pos(2,1)-15/figpos(3) plot_pos(2,2)-.06 plot_pos(2,3)+30/figpos(3) .02];
        maxDNA.Value = (Gates(2,2)-xlims(1))/diff(xlims);
        maxDNA.Visible = 'off';
    end
    manual = uicontrol('style', 'pushbutton');
    manual.Units = 'normalized';
    manual.Position = [.6 .26 .35 .1];
    manual.String = 'Define gate manually';
    manual.Callback = @manualsetup;
    
    approve = uicontrol('style', 'pushbutton');
    approve.Units = 'normalized';
    approve.Position = [.6 .13 .35 .1];
    approve.String = 'Approve';
    approve.Callback = @approveGate;
    
    waitfor(approve, 'backgroundcolor', 'g')
end

%%
    function setGates(src, event, x)
        if x<3
            Gates(1,x) = (diff(ylims)*src.Value)+ylims(1);
        else
            Gates(2,x-2) = (diff(xlims)*src.Value)+xlims(1);
        end
        set(pltgt1, 'XData', [Gates(1,2) max(Gates(1,1), min(p.xLDR))*[1 1] [1 1]*Gates(1,2)])
        if useDNA
            set(pltgt2, 'XData', Gates(2,[1 1 2 2]))
            set(pltgt3, 'XData', Gates(2,[1 1 2 2 1]), 'YData', max(Gates(1,[1 2 2 1 1]),0));
        end
    end

    function manualsetup(src, event)
        minLDR.Visible = 'on';
        maxLDR.Visible = 'on';
        minDNA.Visible = 'on';
        maxDNA.Visible = 'on';
        set(approve, 'backgroundcolor', 'r')
    end

    function approveGate(src, event)
        minLDR.Visible = 'off';
        maxLDR.Visible = 'off';
        minDNA.Visible = 'off';
        maxDNA.Visible = 'off';
        
        set(src, 'backgroundcolor', 'g')
        [LiveCells, DeadCells, Gates, AliveIdx] = DeadCellFilter(LDRtxt, varargin{:}, 'Gates', Gates, 'plotting', false);
    end


end