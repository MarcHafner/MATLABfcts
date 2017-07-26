function [DNAPks, EdUPks] = CCphases(DNA, EdU, plotting)
currfig = gcf;

xDNA = 2.5:.02:8;
xEdU = -.1:.02:5.4;
nsmooth = 5;
EdUshift = 1;
minEdU = 1.5;

if ~exist('plotting','var'), plotting=false; end

% log10 domain and capping
DNA = log10(min(max(DNA, 10^xDNA(2)),10^xDNA(end-1)));
EdU = log10(min(max(EdU, 10^xEdU(2)),10^xEdU(end-1)));

%%
if plotting
    get_newfigure(45654,[50 600 600 650])
    
    % plots with both channels
    subplot(223)
    dscatter(DNA, EdU, 'MSIZE', 20, 'marker', 'o')
    hold on
end

% first pass at the finding peaks by combining DNA and EdU channels

% coarser binning
xDNA2 = xDNA(1:2:end);
xEdU2 = xEdU(1:2:end);
nbins = [length(xDNA2) length(xEdU2)]-1;
bin = NaN(length(DNA),2);
[~,bin(:,2)] = histc(DNA,xDNA2);
[~,bin(:,1)] = histc(EdU,xEdU2);
% counting and smoothing

H = accumarray(bin,1,nbins([2 1]))/length(DNA);
G = smooth1D(H,nsmooth);
F = smooth1D(G',nsmooth)';

% finding peaks
Pk2D = imregionalmax(F);
[x,y] = find(Pk2D);

PksCandidates = [xDNA2(y)' xEdU2(x)' F(Pk2D)];
PksCandidates = sortrows(...
    PksCandidates( (PksCandidates(:,2)>(nsmooth+2)*diff(xEdU2([1 2]))) & ...
    PksCandidates(:,3)>1e-5,:), 3); % filter out the small peaks and one with EdU=0

if plotting
    % plotting the results
    subplot(224)
    imagesc(xDNA2, xEdU2, F)
    hold on
    scatter(PksCandidates(:,1), PksCandidates(:,2), 20+sqrt(PksCandidates(:,3))*100, 'ok')
    set(gca,'ydir','normal')
end

PhasesCandidates = NaN(3,2);

if max(PksCandidates(:,2))-min(PksCandidates(:,2))>EdUshift && ...
        any(PksCandidates(:,2)>(minEdU+.2*EdUshift))
    % there is enough differences on EdU to expect a S phase
    
    % get the S phase peak as the one with max EdU
    PhasesCandidates(2,:) = PksCandidates(argmax(PksCandidates(:,2)),[1 2]);
    
    % find the G1 peak now
    temp = PksCandidates(PksCandidates(:,1)<PhasesCandidates(2,1) & ...
        PksCandidates(:,2)<PhasesCandidates(2,2)-EdUshift, [1 2]);
    if ~isempty(temp) % there is a likely G1 peak
        % assign the G1 peak as the one closest to the S peak on the DNA axis
        PhasesCandidates(1,:) = temp(argmax(temp(:,1)),:);
    end
    % assign the G2 peak as the one closest to the S peak on the DNA axis
    temp = PksCandidates(PksCandidates(:,1)>nanmean(PhasesCandidates(1:2,1)) & ...
        PksCandidates(:,2)<PhasesCandidates(2,2)-EdUshift, [1 2]);
    if ~isempty(temp) % there is a likely G2 peak
        PhasesCandidates(3,:) = temp(argmin(temp(:,1)),:);
    end
    
else
    % most likely no S peak -> take the two highest peak and assign as G1
    % and G2
    PhasesCandidates([1 3],:) = PksCandidates(1:2,[1 2]);
    
end

if plotting
    % plotting the results
    phases = {'G1' 'S' 'G2'};
    for i=1:3
        text(PhasesCandidates(i,1), PhasesCandidates(i,2), phases{i}, ...
            'fontsize', 14, 'fontweight', 'bold', 'color', [.7 .2 .2], ...
            'horizontalalign','center')
    end
end


% now working with each channel sequentially

% work with the DNA content
f = ksdensity(DNA,xDNA);

if plotting
    subplot(221)
    plot(xDNA, f)
    hold on
end

f = ksdensity(DNA(EdU<minEdU+.2*EdUshift),xDNA);
if plotting, plot(xDNA, f, '--'), end

[pks, idx] = findpeaks(f);

% find the DNA peak for G1
DNAPks = idx(sortidx(pks,'descend'));
DNAPks = xDNA(DNAPks(1:min(3,length(DNAPks)))); % take the 3 highest peaks with low EdU
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

if plotting, plot(DNAPks, .1, 'xk'); end


% work with EdU
f = ksdensity(EdU,xEdU);
if plotting
    subplot(222)
    plot(xEdU, f)
    hold on
end
% find the low EdU (G1 and early S)
lE = ( (DNA>DNAPks-1) & (DNA<DNAPks+.1) ) & EdU>2*nsmooth*diff(xEdU(1:2));
f = ksdensity(EdU(lE),xEdU);
if plotting, plot(xEdU, f, '--'), end

[pks, idx] = findpeaks(smooth(f,nsmooth));
EdUPks = idx(sortidx(pks,'descend'));

% take the peak with lowest EdU
EdUPks = min(xEdU(EdUPks([1 2])));


% get the high EdU (S)
hE = (DNA>DNAPks-log10(2)/2) & (DNA<DNAPks+log10(2)*1.5) & EdU>2*nsmooth*diff(xEdU(1:2));
f = ksdensity(EdU(hE),xEdU);
if plotting, plot(xEdU, f, ':'), end

[pks, idx] = findpeaks(smooth(f,nsmooth));
hEdUPks = idx(sortidx(pks,'descend'));
% should be at least 1 above the EdU peak in G1
EdUPks = [EdUPks ...
    xEdU(hEdUPks(find(xEdU(hEdUPks)>(EdUPks+EdUshift),1,'first'))) EdUPks];

% cut off is set half way between the peaks
EdUcutoff = mean(EdUPks);
ylims = [0 EdUPks(2)+(EdUPks(2)-EdUcutoff)];

if plotting
    % plot the location of the peaks
    plot(EdUPks, .1, 'xk');
    plot(EdUcutoff, .1, 'xk');
    xlim(ylims)
end

% get back to the DNA to find the S location
hE = (DNA>DNAPks-log10(2)/2) & (DNA<DNAPks+log10(2)*1.5) & EdU>EdUcutoff;
f = ksdensity(DNA(hE),xDNA);

if plotting
    subplot(221)
    plot(xDNA, f, '-.')
end

[pks, idx] = findpeaks(smooth(f,3*nsmooth));
DNAPks = [DNAPks xDNA(idx(argmax(pks)))];


% get back to the DNA to find the G2
hD = DNA>DNAPks(1)+.3*log10(2) & EdU<EdUcutoff;
f = ksdensity(DNA(hD),xDNA);

if plotting
    subplot(221)
    plot(xDNA, f, ':')
end

[pks, idx] = findpeaks(smooth(f,nsmooth));
hDNAPks = idx(sortidx(pks,'descend'));
% should be around log10(2) above the DNA peak in G1
hDNAPks = xDNA(hDNAPks(xDNA(hDNAPks)>(DNAPks(1)+.5*log10(2))));

if length(hDNAPks)>1
    % more than one candidate
    if ~isnan(PhasesCandidates(3,1))
        % 2D analysis matching G2 peak
        hDNAPks = hDNAPks(argmin(abs(hDNAPks-PhasesCandidates(3,1))));
    elseif ~isnan(PhasesCandidates(2,1)) && any(hDNAPks>PhasesCandidates(2,1))
        % 2D analysis matching S peak
        hDNAPks = min(hDNAPks(hDNAPks>PhasesCandidates(2,1)));
    else
        % take the lowest peak (most likely case in doubt)
        hDNAPks = min(hDNAPks);
    end
end

DNAPks = [DNAPks hDNAPks];

% find the split between G1 and G2
f = ksdensity(DNA,xDNA);
[~, idx] = findpeaks(-smooth(f,5));
DNAcutoff = mean(xDNA(idx(xDNA(idx)>DNAPks(1) & xDNA(idx)<DNAPks(3))));

d1 = DNAcutoff-DNAPks(1);
d2 = DNAPks(3)-DNAcutoff;
xlims = [DNAPks(1)-3*d1 DNAPks(3)+3*d2];
if plotting
    plot(DNAPks, 0, 'xk');
    plot(DNAcutoff, 0, 'xk');
    xlim(xlims)
end


% define areas
G1range = [DNAPks(1)-d1 DNAcutoff];
G2range = [DNAcutoff DNAPks(3)+d2];

if plotting
    % plots with both channels
    subplot(223)
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
    
    
    subplot(224)
    plot(DNAPks, EdUPks, 'xk')
    
    xlim(xlims)
    ylim(ylims)
    pause
    figure(currfig)
end
%%