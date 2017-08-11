function [CCpeaks, CCfrac, Gates, logDNA, logEdU] = CCphases(DNA, EdU, varargin)
currfig = gcf;


p = inputParser;

addParameter(p, 'plotting', false, @islogical)
addParameter(p, 'xDNA', 2.5:.02:8, @isvector)
addParameter(p, 'xEdU', -.2:.02:5.3, @isvector)
addParameter(p, 'nsmooth', 5, @isnumeric)
addParameter(p, 'CCseeds', [], @(x) isempty(x) || (ismatrix(x) && all(size(x)==[3 2])))


parse(p,varargin{:});
p = p.Results;


% determine the spread of the EdU data to calibrate cutoffs
e = -200:1e3;
f = ksdensity(EdU, e);
[~,pk,wdth]=findpeaks(f,'npeaks',1,'widthreference','halfprom','sortstr','descend');
pk = e(ceil(pk));

f2 = ksdensity(EdU(EdU>pk+30), e);
[~,m]=findpeaks(-f2,'npeaks',1,'widthreference','halfprom','sortstr','descend');
m = e(ceil(m));
m = max([m; pk+3*wdth]);


offsetEdU = max(pk-1.5*wdth,1);

% log10 domain and capping
logDNA = log10(min(max(DNA, 10^p.xDNA(3)),10^p.xDNA(end-2)));
logEdU = log10(min(max(EdU-offsetEdU, 10^p.xEdU(3)),10^p.xEdU(end-2)));

% expected maximum EdU value for G1 (in log10)
maxEdU = log10(m-offsetEdU);
% expected minumum EdU value for S (in log10)
minEdU = log10(pk+2*wdth-offsetEdU);
% expected difference between G1/G2 and S (in log10)
EdUshift = max(log10(pk+2*wdth-offsetEdU)-log10(pk-offsetEdU),1);

% % transform the EdU with a sigmoidal
% minEdU = max(quantile(EdU,.02),0)
% logEdU = max(EdU-minEdU,0)./(max(EdU-minEdU,0)+w*2);
% EdUshift = .5; % expected difference between G1/G2 and S (in log10)
% minEdU = .6; % expected minumum EdU value for S (in log10)

%
if p.plotting
    get_newfigure(45654,[5 600 500 750])
    
    % plot original data
    subplot(322)
    idx = randperm(length(EdU),min(length(EdU),1000));
    plot(EdU(idx), logEdU(idx),'.c')
    hold on
    plot(e, max(p.xEdU)*f/max(f), 'k-')
    plot(e, max(p.xEdU)*f2/max(f2), 'k-')
    plot([-200 300 NaN -100 500 ], ...
        [minEdU*[1 1] NaN EdUshift*[1 1]+log10(max(pk-offsetEdU,1)) ], '-r')
    plot(offsetEdU*[1 1], [0 5], '-r')
    plot([-100 m m], [maxEdU*[1 1] 0], '-r')
    ylim(p.xEdU([1 end]))
    xlim([min(EdU) max(pk+5*wdth, 500)])
    
    % plot original data
    subplot(321)
    f = ksdensity(DNA, 0:100:1e6);
    plot(0:100:1e6, f, '-k')
    
    % plots with both channels
    subplot(325)
    dscatter(logDNA, logEdU, 'MSIZE', 20, 'marker', 'o')
    hold on
end

% first pass at the finding peaks by combining DNA and EdU channels

% coarser binning
xDNA2 = p.xDNA(1:2:end);
xEdU2 = p.xEdU(1:2:end);
nbins = [length(xDNA2) length(xEdU2)]-1;
bin = NaN(length(logDNA),2);
[~,bin(:,2)] = histc(logDNA,xDNA2);
[~,bin(:,1)] = histc(logEdU,xEdU2);
% counting and smoothing

H = accumarray(bin,1,nbins([2 1]));
H(H==1) = 0;
H = H/length(logDNA);
G = smooth1D(H,p.nsmooth);
F = smooth1D(G',p.nsmooth)';

% finding peaks
Pk2D = imregionalmax(F);
[x,y] = find(Pk2D);

PksCandidates = [xDNA2(y)' xEdU2(x)' F(Pk2D)];
PksCandidates = sortrows(...
    PksCandidates( (PksCandidates(:,2)>(p.nsmooth+2)*diff(xEdU2([1 2]))) & ...
    PksCandidates(:,3)>1e-5 & PksCandidates(:,3)>max(PksCandidates(:,3)/100),:), 3); % filter out the small peaks and one with EdU=0

if p.plotting
    % plotting the results
    subplot(326)
    imagesc(xDNA2, xEdU2, F)
    hold on
    scatter(PksCandidates(:,1), PksCandidates(:,2), 20+sqrt(PksCandidates(:,3))*100, 'ok')
    set(gca,'ydir','normal')
end

PhasesCandidates = NaN(3,2);

if any(PksCandidates(:,2)-min(PksCandidates(:,2))>EdUshift & ...
        PksCandidates(:,2)>(minEdU+.2*EdUshift))
    % there is enough differences on EdU to expect a S phase
    
    % get the S phase peak as the one with max EdU
    PhasesCandidates(2,:) = PksCandidates(argmax(PksCandidates(:,2)),[1 2]);
    
    % find the G1 peak now
    temp = PksCandidates(PksCandidates(:,1)<PhasesCandidates(2,1)+log10(2)/5 & ...
        PksCandidates(:,1)>PhasesCandidates(2,1)-log10(2) & ...
        PksCandidates(:,2)<PhasesCandidates(2,2)-EdUshift, [1 2]);
    if ~isempty(temp) % there is a likely G1 peak
        % assign the G1 peak as the one closest to the S peak on the DNA axis
        PhasesCandidates(1,:) = temp(argmin(temp(:,1)),:);
    end
    % assign the G2 peak as the one closest to the S peak on the DNA axis
    temp = PksCandidates(PksCandidates(:,1)>nanmean(PhasesCandidates(1:2,1)) & ...
        PksCandidates(:,2)<PhasesCandidates(2,2)-EdUshift, [1 2]);
    if ~isempty(temp) % there is a likely G2 peak
        PhasesCandidates(3,:) = temp(argmin(temp(:,1)),:);
    end
    
else
    % most likely no S peak
    % -> take the two highest peak and assign as G1 and G2
    if ~isempty(PksCandidates)
        PhasesCandidates(1,:) = PksCandidates(1,[1 2]);
        if size(PksCandidates,1)>1 && any((PksCandidates(:,1)-PksCandidates(1,1))>.5*log10(2))
            PhasesCandidates(3,:) = PksCandidates(1+find( ...
                (PksCandidates(2:end,1)-PksCandidates(1,1))>.5*log10(2), 1, 'first'),[1 2]);
        end
    end
end

if p.plotting
    % plotting the results
    phases = {'G1' 'S' 'G2'};
    for i=1:3
        text(PhasesCandidates(i,1), PhasesCandidates(i,2), phases{i}, ...
            'fontsize', 14, 'fontweight', 'bold', 'color', [.7 .2 .2], ...
            'horizontalalign','center')
    end
end

if ~isempty(p.CCseeds)
    PhasesCandidates = p.CCseeds;
    
    if p.plotting
        phases = {'G1' 'S' 'G2'};
        subplot(325)
        for i=1:3
            text(PhasesCandidates(i,1), PhasesCandidates(i,2), phases{i}, ...
                'fontsize', 14, 'fontweight', 'bold', 'color', [0 .8 .8], ...
                'horizontalalign','center')
        end
    
        subplot(326)
        for i=1:3
            text(PhasesCandidates(i,1), PhasesCandidates(i,2), phases{i}, ...
                'fontsize', 14, 'fontweight', 'bold', 'color', [0 .8 .8], ...
                'horizontalalign','center')
        end
    end
end

% now working with each channel sequentially

% work with the DNA content
f = ksdensity(logDNA,p.xDNA);

if p.plotting
    subplot(323)
    plot(p.xDNA, f)
    hold on
end

% take only the cells with low EdU
f = ksdensity(logDNA(logEdU<minEdU+.2*EdUshift & logEdU<maxEdU),p.xDNA);
if p.plotting, plot(p.xDNA, f, '--'), end

[pks, idx] = findpeaks(f, 'sortstr', 'descend');
idx = idx(pks>max(pks/10)); % remove lesser peaks
DNAPks = p.xDNA(idx(1:min(3,length(idx)))); % take the 3 highest peaks with low EdU

% find the DNA peak for G1
if length(DNAPks)>1
    % more than one candidate
    if ~isnan(PhasesCandidates(1,1))
        % 2D analysis matching G1 peak
        DNAPks = DNAPks(argmin(abs(DNAPks-PhasesCandidates(1,1))));
    elseif ~isnan(PhasesCandidates(2,1))
        % 2D analysis matching S peak
        DNAPks = max(DNAPks(DNAPks<PhasesCandidates(2,1)));
    else
        % take the lowest peak (most likely case in doubt)
        DNAPks = min(DNAPks);
    end
end

if p.plotting, plot(DNAPks, .1, 'xk'); end


% work with EdU
f = ksdensity(logEdU,p.xEdU);
if p.plotting
    subplot(324)
    plot(p.xEdU, f)
    hold on
end
% find the low EdU (G1 and early S)
lE = ( (logDNA>DNAPks-1) & (logDNA<DNAPks+.1) ) & ...
    logEdU>2*p.nsmooth*diff(p.xEdU(1:2)) & logEdU<maxEdU;
f = ksdensity(logEdU(lE),p.xEdU);
N = histcounts(logEdU(lE),p.xEdU);
f([true; smooth(N,3)<=1/3]) = 0; % remove single cells
if p.plotting, plot(p.xEdU, f, '--'), end

[~, idx] = findpeaks(smooth(f,p.nsmooth),'sortstr','descend');
% take the highest peak
EdUPks = p.xEdU(idx(1));

%
% get the high EdU (S)
hE = (logDNA>DNAPks-log10(2)/2) & (logDNA<DNAPks+log10(2)*1.5) & logEdU>EdUPks+EdUshift*.8;
if any(hE)
    % some S cells
    f = ksdensity(logEdU(hE),p.xEdU);
    if p.plotting, plot(p.xEdU, f, ':'), end
    
    [pks, idx] = findpeaks(smooth(f,p.nsmooth));
    idx = idx(pks>max(pks/10)); % remove lesser peaks
    pks = pks(pks>max(pks/10));
    hEdUPks = idx(sortidx(pks,'descend'));
    if any(p.xEdU(hEdUPks)>(EdUPks+EdUshift))
        % should be at least 1 above the EdU peak in G1
        EdUPks = [EdUPks ...
            p.xEdU(hEdUPks(find(p.xEdU(hEdUPks)>(EdUPks+EdUshift),1,'first'))) EdUPks];
    else
        % set by default
        EdUPks = [EdUPks (EdUPks+EdUshift) EdUPks];
    end
    
    % cut off is set half way between the peaks
    EdUcutoff = mean(EdUPks([1 2]));
    ylims = [0 min(EdUPks(2)+(EdUPks(2)-EdUcutoff), p.xEdU(end-1))];
    
    if p.plotting
        % plot the location of the peaks
        plot(EdUPks, .1, 'xk');
        plot(EdUcutoff, .1, 'xk');
        xlim(ylims)
    end
    
    % get back to the DNA to find the S location
    hE = (logDNA>DNAPks-log10(2)/2) & (logDNA<DNAPks+log10(2)*1.5) & logEdU>EdUcutoff;
    f = ksdensity(logDNA(hE),p.xDNA);
    
    if p.plotting
        subplot(323)
        plot(p.xDNA, f, '-.')
    end
    
    [pks, idx] = findpeaks(smooth(f,3*p.nsmooth));
    idx = idx(pks>max(pks/10)); % remove lesser peaks
    pks = pks(pks>max(pks/10));
    DNAPks = [DNAPks p.xDNA(idx(argmax(pks)))];
    
else
    % no S cells
    
    % set arbitrary value for high S
    EdUPks = [EdUPks (EdUPks+EdUshift) EdUPks];
    EdUcutoff = mean(EdUPks([1 2]));
    ylims = [0 min(EdUPks(2)+(EdUPks(2)-EdUcutoff), p.xEdU(end-1))];
    
    % set default DNA content
    DNAPks = DNAPks+[0 log10(2)*.5];
end


% get back to the DNA to find the G2
hD = logDNA>DNAPks(1)+.4*log10(2) & logEdU<EdUcutoff;
if any(hD)
    % some G2 cells
    
    f = ksdensity(logDNA(hD),p.xDNA);
    
    if p.plotting
        plot(p.xDNA, f, ':')
    end
    
    [pks, idx] = findpeaks(smooth(f,p.nsmooth), 'sortstr', 'descend');
    idx = idx(pks>max(pks/10)); % remove lesser peaks
    hDNAPks = p.xDNA(idx);
    % should be around log10(2) above the DNA peak in G1
    hDNAPks = hDNAPks(hDNAPks>(DNAPks(1)+.5*log10(2)));
    
    if length(hDNAPks)>1
        % more than one candidate
        if ~isnan(PhasesCandidates(3,1))
            % 2D analysis matching G2 peak
            hDNAPks = hDNAPks(argmin(abs(hDNAPks-PhasesCandidates(3,1))));
        else
            if ~isnan(PhasesCandidates(2,1)) && any(hDNAPks>PhasesCandidates(2,1))
                % 2D analysis matching S peak
                hDNAPks = hDNAPks(hDNAPks>PhasesCandidates(2,1));
            end
            % take the peak closest to a 2-fold(most likely case in doubt)
            hDNAPks = hDNAPks(argmin(abs(hDNAPks-DNAPks(1)-log10(2))));
        end
    end
    
    DNAPks = [DNAPks hDNAPks];
    
    % find the split between G1 and G2
    f = ksdensity(logDNA,p.xDNA);
    [~, idx] = findpeaks(-smooth(f,5));
    DNAcutoff = mean(p.xDNA(idx(p.xDNA(idx)>DNAPks(1) & p.xDNA(idx)<DNAPks(3))));
    if isnan(DNAcutoff)
        DNAcutoff = DNAPks(2);
    end
    
    d1 = DNAcutoff-DNAPks(1);
    d2 = DNAPks(3)-DNAcutoff;
    xlims = [max(DNAPks(1)-3*d1, p.xDNA(2)) min(DNAPks(3)+3*d2, p.xDNA(end-1))];
    if p.plotting
        plot(DNAPks, 0, 'xk');
        plot(DNAcutoff, 0, 'xk');
        xlim(xlims)
    end
else
    % no G2 cells
    
    DNAcutoff = DNAPks(1)+.3*log10(2);
    DNAPks = [DNAPks DNAPks(1)+log10(2)];
    
end

% define areas
G1range = [DNAPks(1)-d1 DNAcutoff];
G2range = [DNAcutoff DNAPks(3)+d2];
Gates = struct('G1', [G1range' [p.xEdU(1) EdUcutoff]'], ...
    'S', [ [G1range(1) G2range(2)]', [EdUcutoff EdUPks(2)+EdUcutoff*.5]'], ...
    'G2', [G2range' [p.xEdU(1) EdUcutoff]'] );

CCfrac = [mean(logDNA>=G1range(1) & logDNA<G1range(2) & logEdU<EdUcutoff) ...
    mean(logDNA>=G1range(1) & logDNA<G2range(2) & logEdU>=EdUcutoff) ...
    mean(logDNA>=G2range(1) & logDNA<G2range(2) & logEdU<EdUcutoff) ...
    mean(logDNA<G1range(1) | logDNA>G2range(2) )];

CCpeaks = [DNAPks' EdUPks'];

if p.plotting
    % plots with both channels
    subplot(325)
    plot(G1range([1 1 2 2]), [0 EdUcutoff EdUcutoff 0], '-r')
    plot(G2range([1 1 2 2]), [0 EdUcutoff EdUcutoff 0], '-r')
    plot([G1range([1 1]) G2range([2 2])], ...
        [EdUPks(2)+EdUcutoff*.5 EdUcutoff EdUcutoff EdUPks(2)+EdUcutoff*.5], '-r')
    
    phases = {'G1' 'S' 'G2'};
    for i=1:3
        text(DNAPks(i), EdUPks(i), phases{i}, ...
            'fontsize', 14, 'fontweight', 'bold', 'color', 'k', ...
            'horizontalalign','center')
    end
    xlim(xlims)
    ylim(ylims)
    
    
    subplot(326)
    plot(DNAPks, EdUPks, 'xk')
    
    xlim(xlims)
    ylim(ylims)
end
%%
if p.plotting
    pause
    figure(currfig)
end