function [CCfrac, pH3CellIdentity, pH3cutoff] = pH3Filter(pH3, CellIdentity, varargin)

assert(all(size(pH3)==size(CellIdentity)))

p = inputParser;

addParameter(p, 'plotting', false, @islogical)
addParameter(p, 'interactive', false, @islogical)
addParameter(p, 'xpH3', 2.5:.02:8, @isvector)
addParameter(p, 'pH3cutoff', [], @isscalar)
% addParameter(p, 'Gates', NaN(2), @(x) ismatrix(x) & all(size(x)==2) & all(~isnan(x(:,1))))
addParameter(p, 'savefigure', '', @ischar)

parse(p,varargin{:});
p = p.Results;

if p.interactive, p.plotting = true; end
if ~isempty(p.savefigure), p.plotting = true; end


%%


logpH3 = log10(min(max(pH3, 10^p.xpH3(3)),10^p.xpH3(end-2)));

% only take the cells in G1 or G2
f = ksdensity(logpH3(CellIdentity==1 | CellIdentity==3), p.xpH3, 'width', 2.5*diff(p.xpH3(1:2)));

if ~isempty(p.pH3cutoff)
    % already defined gate
    pH3cutoff = p.pH3cutoff;
else
    % determine the spread of the pH3 data to define cutoff
    [~, pk, pH3wdth] = findpeaks(f,'npeaks',1,'widthreference','halfprom','sortstr','descend');
    
    [~, minpk] = findpeaks(-f(pk:end),'npeaks',1); minpk=minpk+pk-1;
    pH3cutoff = p.xpH3(ceil(max(min(minpk, pk+6*pH3wdth), pk+2*pH3wdth)));
end


if p.plotting
    % plotting results
    currfig = gcf;
    
    plot_pos = [
    .07 .3 .4 .67;
    .6 .55 .3 .4
    .6 .08 .3 .4];

    get_newfigure(45679,[505 100 550 300])
    xlims = quantile(logpH3, [1e-3 .999]);
    
    % plot original data
    get_newaxes(plot_pos(1,:),1)
    plot(p.xpH3, log10(f+max(f)/100)-log10(max(f)/100))
    
    fall = ksdensity(logpH3, p.xpH3, 'width', 2.5*diff(p.xpH3(1:2)));
    plot(p.xpH3, log10(fall+max(fall)/100)-log10(max(fall)/100), '--')
    
    plot([1 1]*pH3cutoff, ...
        [0 .5]*log10(max(f)), '-', 'color', [.6 .9 .1]);
    pltgt1 = plot([1 1]*pH3cutoff, ...
        [0 .5]*log10(max(f)), 'r-');
    
    xlim(xlims)
    ylim([0 log10(max(f))-log10(max(f)/100)+.1])
    set(gca,'xtick',[],'ytick',[])
    
    pieax = get_newaxes(plot_pos(2,:));   
    set(gca,'xtick',[],'ytick',[],'visible','off') 
    axis square
    
    pieax2 = get_newaxes(plot_pos(3,:));   
    set(gca,'xtick',[],'ytick',[],'visible','off') 
    axis square
end


% finalize the results
EvalMphase();



if p.interactive
    figpos = get(gcf,'position');
    
    minpH3 = uicontrol('style', 'slider', 'callback', {@setCutoff,1});
    minpH3.Units = 'normalized';
    minpH3.InnerPosition = [plot_pos(1,1)-15/figpos(3) plot_pos(1,2)-.04 plot_pos(1,3)+30/figpos(3) .03];
    minpH3.Value = max(0,(pH3cutoff-xlims(1))/diff(xlims));
    
    
    
    approve = uicontrol('style', 'pushbutton');
    approve.Units = 'normalized';
    approve.Position = [.1 .03 .3 .13];
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
    function setCutoff(src, event, x)
        pH3cutoff = (diff(xlims)*src.Value)+xlims(1);        
        set(pltgt1, 'XData', pH3cutoff*[1 1])
        EvalMphase();
    end


    function Midx = EvalMphase()
        Midx = logpH3>=pH3cutoff;
        pH3CellIdentity = CellIdentity;
        pH3CellIdentity(Midx & ismember(CellIdentity,[1 3])) = 4;        
        
        for id = 1:5
            CCfrac(id) = mean(pH3CellIdentity==mod(id,5));
        end
        
        if p.plotting
            set(gcf,'currentaxes', pieax)
            cla
            ptxt = pie(mean([Midx ~Midx]) ,{sprintf('M-phase cells (%.1f%%)', ...
                100*mean(Midx)), 'Other'});
            set(ptxt(2),'fontsize',12, 'fontweight','bold')
            
            set(gcf,'currentaxes', pieax2)
            cla
            pie(CCfrac, {'G1' 'S' 'G2' 'M' 'other'})
        end
    end

    function approveGate(src, event)
        minpH3.Visible = 'off';
        
        set(src, 'backgroundcolor', 'g')
    end


end