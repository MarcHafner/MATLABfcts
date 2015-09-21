function MOMP_data = Evaluate_MOMPTime(Raw_data, plotting, verbatim)
% MOMP_data = Evaluate_MOMPTime(Raw_data, plotting, verbatim)
%
%   MOMP_data.idx_momp.Edge (and MOMP_data.idx_momp.MOMP)
%                   = -1 for too short,
%                    Inf for survivors,
%                 frame# for MOMP at Time(frame#)
%                     -2 for unclassified
%
%   MOMP_data.idx_momp.lostidx = idx for tracking lost.
%
%   MOMP_data.Traj.Edge = edge (smoothed)
%   MOMP_data.Traj.Area = area (smoothed)
%   MOMP_data.Traj.Prob = edge/area scoring
%

fprintf('\n---- Function Evaluate_MOMPTime for %s -----\n\n', Raw_data.tag)

if ~exist('plotting','var') || isempty(plotting)
    plotting = 0;
end
if ~exist('verbatim','var')
    verbatim = false;
end


windowSize = round(45/Raw_data.Timestep);         % windows size for filtering
division_time = round(250/Raw_data.Timestep);     % assumed division time
edge_cutoff = 100;      % cutoff for the edge criterion
MOMPchannel_cutoff = 70;          % cutoff for the MOMP channel
margin_delta = .1;      % hyteresis for the scoring of the derivate
p_cutoff_high = 1.1;    % lower cutoff for the score to be classified as MOMP
p_cutoff_low = .5;      % upper cutoff for the score to be classified as survival
MOMP_cutoff = Raw_data.Ntimepoints-windowSize;
Delta_Death_MOMP = 4;   % number of frame for MOMP to occur before cell change shape

MOMP_data = Raw_data;
MOMP_data.parameters = struct(...
    'MOMPwindowSize', windowSize, ...
    'division_time', division_time, ...
    'edge_cutoff', edge_cutoff, ...
    'MOMPchannel_cutoff', MOMPchannel_cutoff, ...
    'margin_delta', margin_delta, ...
    'p_cutoff_high', p_cutoff_high, ...
    'p_cutoff_low', p_cutoff_low, ...
    'MOMP_cutoff', MOMP_cutoff, ...
    'Delta_Death_MOMP', Delta_Death_MOMP);

if strcmpi(MOMP_data.RFPchannel, 'MOMP')
    %     assert(all(cell2mat( cellfun2(@(x) any(~isnan(x)), {MOMP_data.rawdata.momp}))), ...
    %         'some trajectories have NaN in the momp variable whereas MOMP was specified to be used')
    useMOMPchannel = true;
else
    useMOMPchannel = false;
end

idx_momp = struct('MOMPchannel',{},'Edge',{},'track_lost',{});
Traj = struct('Prob',{},'Edge',{},'MOMPchannel',{},'Area',{});


WellLabel = MOMP_data.WellLabel;
Wells = MOMP_data.Wells(:,3);
rawdata = MOMP_data.rawdata;
parameters = MOMP_data.parameters;

 parfor (iW=1:MOMP_data.Nwells, (plotting==0)*100)
% for iW=1:MOMP_data.Nwells
    fprintf('\n*** Analyzing the well %s (%s) ***\n', Wells{iW},WellLabel{iW})

    [idx_momp(iW), Traj(iW)] = evalMOMP(rawdata(iW), ...
        plotting, parameters, useMOMPchannel, verbatim);

    fprintf('\tSurvivors: %i, MOMP: %i, unclassified: %i, too short: %i\n', ...
        sum(isinf(idx_momp(iW).Edge)), ...
        sum((idx_momp(iW).Edge>0) & (idx_momp(iW).Edge<Inf)), ...
        sum(idx_momp(iW).Edge==-2), sum(idx_momp(iW).Edge==-1));
end

MOMP_data.idx_momp = idx_momp;
MOMP_data.Traj = Traj;

end


function [MOMPtime, traj] = evalMOMP(rawdata, plotting, para, useMOMPchannel, verbatim)

windowSize = para.MOMPwindowSize;
division_time = para.division_time;
edge_cutoff = para.edge_cutoff;
MOMPchannel_cutoff = para.MOMPchannel_cutoff;
margin_delta = para.margin_delta;
p_cutoff_high = para.p_cutoff_high;
p_cutoff_low = para.p_cutoff_low;
MOMP_cutoff = para.MOMP_cutoff;
Delta_Death_MOMP = para.Delta_Death_MOMP;


n_cell = size(rawdata.fret, 2);

MOMPtime_Edge = NaN(n_cell, 1);
MOMPtime_track_lost = NaN(n_cell, 1);
traj_Edge = NaN(MOMP_cutoff, n_cell);
traj_Area = NaN(MOMP_cutoff, n_cell);
traj_Prob = NaN(MOMP_cutoff, n_cell);
MOMPtime_MOMPchannel = NaN(n_cell, 1);
traj_MOMPchannel = NaN(MOMP_cutoff, n_cell);

temp_plot = 'y';

edge = rawdata.edge;
area = rawdata.area;
momp = rawdata.momp;
cfp = rawdata.cfp;
fret = rawdata.fret;

fprintf('%3i cells to analyze:                      |', n_cell)
%  parfor (q = max(min(plotting*(length(plotting)>1)),1) : n_cell, plotting)
for q = max(min(plotting*(length(plotting)>1)),1) : n_cell

    if ~verbatim && mod(q,200)==1
        fprintf('\n%3i: ', q-1)
    elseif ~verbatim && mod(q,5)==1
        fprintf('.');
    end

    e = edge(:, q);
    a = area(:, q);
    if useMOMPchannel
        m = momp(:, q);end

    idx_last = find(~isnan(e), 1, 'last');

    if isempty(idx_last) || (~isempty(idx_last) && idx_last < windowSize)
        if verbatim, fprintf('\t%d too short\n', q);end
        MOMPtime_Edge(q) = -1;
        if useMOMPchannel
            MOMPtime_MOMPchannel(q) = -1;end
        continue
    elseif ~isempty(idx_last)
        MOMPtime_track_lost(q) = idx_last;
    else
        error('unexpected case');
    end

    idx_last = min(idx_last, MOMP_cutoff+floor(windowSize/2));


    E = filter_smooth(e, idx_last, windowSize);
    A = filter_smooth(a, idx_last, windowSize);
    traj_Edge(:, q) = [E(1:min(idx_last,MOMP_cutoff));
        NaN(MOMP_cutoff-min(idx_last,MOMP_cutoff),1)];
    traj_Area(:, q) = [A(1:min(idx_last,MOMP_cutoff));
        NaN(MOMP_cutoff-min(idx_last,MOMP_cutoff),1)];
    if useMOMPchannel
        M = filter_smooth(m,idx_last, windowSize);
        traj_MOMPchannel(:, q) = [M(1:min(idx_last,MOMP_cutoff));
            NaN(MOMP_cutoff-min(idx_last,MOMP_cutoff),1)];
    end

    %%

    sliding_edgecutoff = edge_cutoff+nanmedian(E(5+[1:windowSize]));
    p = zeros(4,length(e)); % based on edge cutoff, edge difference, area difference and tracking
    p2 = zeros(1,length(e));

    noise_idx = 4+find(e(5:end)>sliding_edgecutoff,1,'first');
    Enoise = eval_Trajnoise(e, E, noise_idx, windowSize);
    Anoise = eval_Trajnoise(a, A, noise_idx, windowSize);


    if useMOMPchannel
        sliding_edgecutoff = MOMPchannel_cutoff+nanmedian(M(5+[1:windowSize]));
        Mnoise_idx = 4+find(m(5:end)>sliding_edgecutoff,1,'first');
        p_momp = zeros(3,length(e)); % based on MOMP cutoff, MOMP difference and tracking
        p2_momp = zeros(1,length(e));
        temp_Mcutoff = MOMPchannel_cutoff+nanmean(M(1:windowSize));
        Mnoise = eval_Trajnoise(m, M, Mnoise_idx, windowSize);
    end

    % score based on the derivate of the edge
    [p_deltaE, deltaE] = derivate_score(E, idx_last, windowSize, edge_cutoff, MOMP_cutoff, margin_delta);
    p(4,1:length(p_deltaE)) = .4*p_deltaE;

    if useMOMPchannel
        % score based on the derivate of the MOMP channel
        p_deltaM = derivate_score(M, idx_last, windowSize, edge_cutoff, MOMP_cutoff, margin_delta);
        p_momp(3,1:length(p_deltaM)) = .4*p_deltaM;
    end

    for i=windowSize:idx_last
        % different windows
        pre_2w = [max(windowSize-2,i-2*windowSize):(i-1)];
        pre_w = [max(windowSize-2,i-windowSize+1):(i-1)];
        post_w = [i:min([max(length(E)-5,i+windowSize), i+2*windowSize, length(E)])];
        death_w = [i:min([i+division_time, max(i+5,length(E)-5), length(E)-1])];

        if sliding_edgecutoff > nanmean(E(pre_w))+Enoise/2
            sliding_edgecutoff = edge_cutoff+nanmedian(E(pre_2w));
            % if verbatim, fprintf('\t  edgecutoff update at %i to %.2f\n', i, sliding_edgecutoff);end
        end

        % criteria based on edge: cutoff and ranksum test
        prob = cutoff_dist_prob(E, sliding_edgecutoff, edge_cutoff, Enoise, pre_w, pre_2w, post_w, death_w);
        p(1:2,i) = [.8 .4]'.*prob;
        % criteria based on area
        %                 p(3,i) = .3* max(0,(.05 -cdf('normal',max(A(death_w)), mean(A(pre_2w-2)), Anoise+std(A(pre_2w-2))))/.05) ...
        %                     -.2* max(0,(cdf('normal',max(A(death_w)), mean(A(i:max(post_w-2))), Anoise+std(A(i:max(post_w-2)))) -.8)/.2);
        p(3,i) = .3* max(0,(.05 -precomputed_normal_cdf(max(A(death_w)), mean(A(pre_2w-2)), Anoise+std(A(pre_2w-2))))/.05) ...
            -.2* max(0,(precomputed_normal_cdf(max(A(death_w)), mean(A(i:max(post_w-2))), Anoise+std(A(i:max(post_w-2)))) -.8)/.2);

        % contribution for end of tracking
        p(4,i) = p(4,i)* (1+ .5*max(0, min(1, ...
            (windowSize+1 -(idx_last-i-5))/(windowSize+1)))*(idx_last<MOMP_cutoff));

        p2(i) = sum(p(:,i))*(E(i)>sliding_edgecutoff-Enoise);

        if useMOMPchannel
            prob = cutoff_dist_prob(M, temp_Mcutoff, 0, Mnoise, pre_w, pre_2w, post_w, death_w);
            p_momp(1:2,i) = [.8 .4]'.*prob;

            p_momp(3,i) = p_momp(3,i)* (1+ .5*max(0, min(1, ...
                (windowSize+1 -(idx_last-i-5))/(windowSize+1)))*(idx_last<MOMP_cutoff));

            p_momp(:,i) = p_momp(:,i)*1.9/1.6; % rescaling because not using area
            p2_momp(i) = sum(p_momp(:,i))*(M(i)>temp_Mcutoff-Mnoise);
        end
    end

    MOMPtime_Edge(q) = DetermineMOMP(p, p2, Delta_Death_MOMP, MOMP_cutoff, p_cutoff_high, p_cutoff_low);
    traj_Prob(:,q) = [p2(1:min(length(p2),MOMP_cutoff));
        NaN(MOMP_cutoff-min(length(p2),MOMP_cutoff),1)];


    if verbatim
        if isinf(MOMPtime_Edge(q))
            fprintf('\t%d MOMP at %i\n', q, MOMPtime_Edge(q))
        elseif MOMPtime_Edge(q)==0
            fprintf('\t%d survive\n', q)
        else
            fprintf('\t%d unclassified\n', q);
        end
    end
    if useMOMPchannel
        MOMPtime_MOMPchannel(q) = DetermineMOMP(p_momp, p2_momp, 0, MOMP_cutoff, p_cutoff_high, p_cutoff_low);
    end
    %%

    if any(plotting~=0) && temp_plot=='y' && ...
            ((length(plotting)>1 && ismember(q,abs(plotting)) || ...
            (length(plotting)==1 && q<abs(plotting))))
        %%
        figure(993);clf
        if ~useMOMPchannel
            m = NaN*e;
            M = NaN*e;
            tM = NaN*traj_Edge(:,q);
        end

        subplot(3,2,[1 3])
        hold on
        plot(m,'r','linewidth',1)
        plot(a,'g','linewidth',1)
        plot(A,'g','linewidth',2)
        plot(e,'b','linewidth',1)
        plot(E,'b','linewidth',2)
        plot(150+deltaE*1000,'c','linewidth',2)

        if MOMPtime_Edge(q)>0
            plot(MOMPtime_Edge(q)*[1 1], [0 max([e;a])], '.k-', 'linewidth', 2)
            title('MOMP detected')
        elseif MOMPtime_Edge(q)==0
            title('survival')
        elseif MOMPtime_Edge(q)==-1
            title('Too short')
        elseif MOMPtime_Edge(q)==-2
            title('Uncertain classification')
        end

        if useMOMPchannel && MOMPtime_Edge(q)>0
            plot(MOMPtime_Edge(q)*[1 1], [0 max([e;a])], '.r-','markersize',20)
        end
        xlim([0 length(e)])
        ylim([0 nanmax([e;a])])

        subplot(3,2,5)
        hold on
        plot(p(1,:), 'b')
        plot(p(2,:), 'b:')
        plot(p(3,:), 'g')
        plot(p(4,:), 'm')
        plot(p2, 'k:', 'linewidth', 2)
        if MOMPtime_Edge(q)>0
            plot(MOMPtime_Edge(q)*[1 1], [0 2], '.k-', 'linewidth', 2)
            %                     title(['Max p = ' num2str(p2(MOMPtime_Edge(q)))])
        else
            %                     title(['Max p = ' num2str(max(sum(p)))])
        end
        if useMOMPchannel && MOMPtime_MOMPchannel(q)>0
            plot(MOMPtime_MOMPchannel(q)*[1 1], [0 2], '.r-','markersize',20)
        end

        xlim([0 length(e)])
        ylim([-.25 2])


        subplot(122)
        hold on
        plot(cfp(:,q)./fret(:,q),'g','linewidth',1)
        plot(cfp(:,q),'c','linewidth',1)
        plot(100+fret(:,q)*200,'b','linewidth',2)
        if MOMPtime_Edge(q)~=0
            plot(MOMPtime_Edge(q)*[1 1], [0 max(cfp(:,q))], '.k-')
        end
        if useMOMPchannel && MOMPtime_MOMPchannel(q)>0
            plot(MOMPtime_MOMPchannel(q)*[1 1], [0 max(cfp(:,q))], '.r-','markersize',20)
        end

        xlim([0 length(e)])
        %%

        if plotting>0
            temp_plot = input('    Do you want more? y/n [y]: ', 's');
            if isempty(temp_plot)
                temp_plot = 'y';
            end
        else
            pause(.3)
        end
    end


end
MOMPtime.MOMPchannel = MOMPtime_MOMPchannel;
MOMPtime.Edge = MOMPtime_Edge;
MOMPtime.track_lost = MOMPtime_track_lost;

traj.Prob = traj_Prob;
traj.Edge = traj_Edge;
traj.MOMPchannel = traj_MOMPchannel;
traj.Area = traj_Area;

fprintf('\n');
end




function x = filter_smooth(x, idx_last, windowSize)
if idx_last>(3*windowSize)
    x(1:idx_last) = filtfilt(ones(1,windowSize)/windowSize,1,x(1:idx_last));
else
    x(1:idx_last) = smooth(x(1:idx_last), windowSize, 'rlowess');
end
end

function noise = eval_Trajnoise(x, smoothx, noise_idx, windowSize)
noise = nanstd(x(5:min(5+windowSize,noise_idx)) - ...
    smoothx(5:min(5+windowSize,noise_idx)));
noise = nanmin(noise, 100);
end


function [p, delta] = derivate_score(x, idx_last, windowSize, edge_cutoff, MOMP_cutoff, margin_delta)
delta = max(-.15,[0; diff((x-nanmean(x(1:windowSize)))/edge_cutoff)]);
p = zeros(1,min(idx_last, MOMP_cutoff-floor(windowSize/2)));
for i=2:length(p)
    if delta(i)>margin_delta
        p(i) = min(1,max(0,p(i-1))+.05+min(.15,(delta(i)-margin_delta)/margin_delta));
    elseif delta(i)<-margin_delta
        p(i) = max(-.3, min(0,p(i-1))-.03);
    else
        p(i) = 0;
    end
end
end


function prob = cutoff_dist_prob(x, sliding_cutoff, cutoff, xnoise, pre_w, pre_2w, post_w, death_w)
prob = max(0, min(1, ...
    sum(x(post_w)>sliding_cutoff+xnoise)/max(length(post_w),5) ...
    - .5*sum(x(pre_w)>sliding_cutoff+.2*cutoff+xnoise)/max(length(pre_w),5) ...
    - .25*sum(x(pre_w)>(sliding_cutoff+1.5*cutoff+xnoise))/max(length(pre_w),5) ...
    - sum(x(death_w)<sliding_cutoff)/max(length(death_w),5)));

% prob(2,1) = max(0,(cdf('normal',min(x(death_w)), mean(x(pre_2w-2)), xnoise+std(x(pre_2w-2))) -.9)/.1);
prob(2,1) = max(0,(precomputed_normal_cdf(min(x(death_w)), mean(x(pre_2w-2)), xnoise+std(x(pre_2w-2))) -.9)/.1);
end


function MOMPtime = DetermineMOMP(prob_MOMP, prob2_MOMP, delta_D_MOMP, MOMP_cutoff, p_cutoff_high, p_cutoff_low)
if any(prob2_MOMP(1:MOMP_cutoff)>p_cutoff_high)
    temp = find(prob2_MOMP(1:MOMP_cutoff)>p_cutoff_high,1,'last');
    MOMP_range = find(prob2_MOMP(1:temp)>p_cutoff_high, 1, 'first'):temp;

    MOMPtime = find( (prob2_MOMP(MOMP_range) > max(prob2_MOMP(MOMP_range))*.95) & ...
        (prob2_MOMP(MOMP_range)>p_cutoff_high) , 1, 'last');
    MOMPtime = MOMPtime+MOMP_range(1)-1 -delta_D_MOMP;

elseif all(sum(prob_MOMP(:,1:MOMP_cutoff))<p_cutoff_low)
    MOMPtime = Inf;
else
    MOMPtime = -2;
end
end
