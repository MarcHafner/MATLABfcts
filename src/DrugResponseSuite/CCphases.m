function [nodes, boundaries] = CCphases(DNA, EdU)
currfig = gcf;

xDNA = 3:.02:8;
xEdU = -.06:.02:5;

% should be in the log10 domain;
DNA = log10(min(max(DNA, 10^xDNA(1)),10^xDNA(end)));
EdU = log10(min(max(EdU, 10^xEdU(1)),10^xEdU(end)));

%%
get_newfigure(45654,[50 600 600 650])


f = ksdensity(DNA,xDNA);

subplot(221)
plot(xDNA, f)
hold on

[pks, idx] = findpeaks(f);
PhasePks = idx(sortidx(pks,'descend'));


PhasePks = sort([PhasePks(1) ...
    PhasePks(1+find((abs(xDNA(PhasePks(2:min(end,10)))-xDNA(PhasePks(1)))<log10(2.8)) & ...
    (abs(xDNA(PhasePks(2:min(end,10)))-xDNA(PhasePks(1)))>log10(1.4)),...
    1,'first'))]);


[~, idx] = findpeaks(-smooth(f,5));
DNAcutoff = mean(xDNA(idx(idx>PhasePks(1) & idx<PhasePks(2))));

PhasePks = xDNA(PhasePks);
d1 = DNAcutoff-PhasePks(1);
d2 = PhasePks(2)-DNAcutoff;
plot(PhasePks, 0, 'xk');
plot(DNAcutoff, 0, 'xk');
xlim([PhasePks(1)-3*(DNAcutoff-PhasePks(1)) PhasePks(2)+3*(PhasePks(2)-DNAcutoff)])


% get the split for the EdU

f = ksdensity(EdU,xEdU);
subplot(222)
plot(xEdU, f)
hold on

% low EdU (G1 and early S)
lE = ( (DNA>PhasePks(1)-d1) & (DNA<PhasePks(1)+.3*d1) ) & EdU>.13;
f = ksdensity(EdU(lE),xEdU);
plot(xEdU, f, '--')

[pks, idx] = findpeaks(smooth(f,5));
EdUPks = idx(sortidx(pks,'descend'));

% take the peak with lowest EdU
EdUPks = min(xEdU(EdUPks([1 2])));

% high EdU (S)
hE = (DNA>PhasePks(1)+.4*d1) & (DNA<PhasePks(2)-.4*d2) & EdU>.13;
f = ksdensity(EdU(hE),xEdU);
plot(xEdU, f, ':')

[pks, idx] = findpeaks(smooth(f,5));
hEdUPks = idx(sortidx(pks,'descend'));
EdUPks = [EdUPks ...
    xEdU(hEdUPks(find(xEdU(hEdUPks)>(EdUPks+1),1,'first')))];

EdUcutoff = mean(EdUPks);

plot(EdUPks, 0, 'xk');
plot(EdUcutoff, 0, 'xk');
xlim([0 EdUPks(2)+EdUcutoff])

% define areas
G1range = [PhasePks(1)-(DNAcutoff-PhasePks(1)) DNAcutoff];
G2range = [DNAcutoff PhasePks(2)+(PhasePks(2)-DNAcutoff)];

subplot(223)
dscatter(DNA, EdU, 'MSIZE', 20, 'marker', 'o')
hold on
plot(G1range([1 1 2 2]), [0 EdUcutoff EdUcutoff 0], '-r')
plot(G2range([1 1 2 2]), [0 EdUcutoff EdUcutoff 0], '-r')
plot([G1range([1 1]) G2range([2 2])], ...
    [EdUPks(2)+EdUcutoff*.5 EdUcutoff EdUcutoff EdUPks(2)+EdUcutoff*.5], '-r')
plot([PhasePks(1) DNAcutoff PhasePks(2)], EdUPks([1 2 1]), 'xk')

xlim([PhasePks(1)-3*(DNAcutoff-PhasePks(1)) PhasePks(2)+3*(PhasePks(2)-DNAcutoff)])
ylim([0 EdUPks(2)+EdUcutoff])


%
xDNA2 = xDNA(1:2:end);
xEdU2 = xEdU(1:2:end);
nbins = [length(xDNA2) length(xEdU2)]-1;
bin = NaN(length(DNA),2);
[~,bin(:,2)] = histc(DNA,xDNA2);
[~,bin(:,1)] = histc(EdU,xEdU2);
H = accumarray(bin,1,nbins([2 1]));
G = smooth1D(H,4);
F = smooth1D(G',4)';

Pk2D = imregionalmax(F);
[x,y] = find(Pk2D');

PksCandidates = [xDNA2(x), xEdU2(y)];

subplot(224)    
imagesc(xDNA2, xEdU2, F)
hold on
plot(xDNA2(x), xEdU2(y), 'xk')
xlim([PhasePks(1)-3*(DNAcutoff-PhasePks(1)) PhasePks(2)+3*(PhasePks(2)-DNAcutoff)])
ylim([0 EdUPks(2)+EdUcutoff])
set(gca,'ydir','normal')
%%
pause
figure(currfig)