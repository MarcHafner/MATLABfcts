function Process_data = Process_trajectories(MOMP_data, plotting)
% Process_data = Process_trajectories(MOMP_data, plotting)
%
%
%   Using the edge as a marker of cell death.
%   Cutoff at maxT (10 hours)
%
% case 0      good tracking: survival until 10h or MOMP before end of tracking
% case 1      Bad tracking
% case 2      Tracking lost before minT(5h), no MOMP
% case 12     Tracking lost between minT(5h) and maxT(10h), no MOMP
% case 3      Undetermined fate
% case 4      too short, less than minlength (10) frames
% case 5      Unexpected case
%
fprintf('\n---- Function Process_trajectories for %s -----\n\n', MOMP_data.tag)

if ~exist('plotting','var')
    plotting = false;
end

Process_data = MOMP_data;
Process_data.fullT = (1:Process_data.Ntimepoints)*Process_data.Timestep;
Process_data.parameters.smoothing = max(Process_data.parameters.MOMPwindowSize-2,5);

windowSize = round(Process_data.parameters.MOMPwindowSize+10/Process_data.Timestep);
minlength = round(50/Process_data.Timestep); % minimum number of frame for considering a good tracking
maxT = 600; % number of minutes for full tracking
minT = 300; % for a full tracking except for cell death
TafterMOMP = 600; % Time after MOMP for trajectories
SubstractedFrames = 4;

Process_data.parameters.ProcesswindowSize = windowSize;
Process_data.parameters.minlength = minlength;
Process_data.parameters.maxT = maxT;
Process_data.parameters.minT = minT;
Process_data.parameters.TafterMOMP = TafterMOMP;
Process_data.parameters.SubstractedFrames = SubstractedFrames;


% controls
DrugIdx = ~ismember(Process_data.Drugs, {'FlipL' 'FlipS' 'BCL-2' 'BCL-XL'});
ctrlidx = find(all(Process_data.Doses(:,DrugIdx)==0,2));
figure(993);clf
ctrl = NaN(Process_data.Ntimepoints, length(ctrlidx));
for i=1:length(ctrlidx)
    ctrl(:,i) = smooth(smooth(trimmean(Process_data.rawdata(ctrlidx(i)).fret,50,2),...
        Process_data.parameters.smoothing,'rlowess'),2*Process_data.parameters.smoothing,'rlowess');

    subplot(floor(sqrt(length(ctrlidx))), ...
        ceil(length(ctrlidx)/floor(sqrt(length(ctrlidx)))), i)
    plot(Process_data.fullT,Process_data.rawdata(ctrlidx(i)).fret)
    hold on
    plot(Process_data.fullT,ctrl(:,i),'-k','linewidth',2)
    hold off
end
ctrl = median(ctrl,2);
Process_data.ctrl = ctrl;


% Determine the cutoff (in minutes) to consider survival
cutoff = find(Process_data.fullT>=maxT, 1, 'first');
cutoff2 = cutoff +windowSize +2;
Process_data.T = Process_data.fullT(1:cutoff);
T2 = Process_data.fullT(1:cutoff2); % slightly longer T for smoothing

warning('********************')
    warning('iW not complete')
for iW=1:Process_data.Nwells
    fprintf('\n** Analyzing the well %s (%s) **\n', Process_data.Wells{iW,3}, ...
        Process_data.WellLabel{iW})

    % application of the survival cutoff
    FRET = Process_data.rawdata(iW).fret(1:cutoff2,:) ...
        -Process_data.ctrl(1:cutoff2)*ones(1,Process_data.Ntrajectories(iW));
    Process_data.idx_momp(iW).Edge(Process_data.idx_momp(iW).Edge>=cutoff) = Inf;   % surviving cells
    Process_data.idx_momp(iW).track_lost(Process_data.idx_momp(iW).track_lost>=cutoff) = NaN;    % surviving cells

    idx_momp = Process_data.idx_momp(iW).Edge;

    % remove the bad tracks
    goodtrack = false(1,Process_data.Ntrajectories(iW));
    ErrorType = zeros(1,Process_data.Ntrajectories(iW));

    % categorize the 'bad tracks'
    Process_data.Legends.ErrorType = containers.Map(...
        [0 1 2 12 3 4 5],{'Good tracking', 'Bad tracking', ...
        'Tracking lost before 5h, no MOMP',...
        sprintf('Tracking lost between %.1fh and %.1fh, no MOMP', minT/60, maxT/60), ...
        'Undetermined fate',sprintf('too short (less than %i frames)',minlength),...
        'Unexpected case'});

    for j=1:Process_data.Ntrajectories(iW)

        if Process_data.idx_momp(iW).Edge(j)==-1   % -1 = bad tracking
            ErrorType(j) = 1;
        elseif Process_data.idx_momp(iW).Edge(j)==-2   % -2 = undetermined fate
            ErrorType(j) = 3;
        elseif isinf(Process_data.idx_momp(iW).Edge(j))
            % Inf = survive => need to be tracked until the end
            goodtrack(j) = isnan(Process_data.idx_momp(iW).track_lost(j));
            if ~goodtrack(j)
                if T2(Process_data.idx_momp(iW).track_lost(j))>minT
                    ErrorType(j) = 12; % tracked to the minimum time
                else
                    ErrorType(j) = 2; % tracking too short
                end
            end
        elseif Process_data.idx_momp(iW).Edge(j)>0 && Process_data.idx_momp(iW).Edge(j)<Inf && ...
                (Process_data.idx_momp(iW).track_lost(j)>=minlength || ...
                isnan(Process_data.idx_momp(iW).track_lost(j)) )
            % momp recorded => need to be tracked until momp at least
            if (Process_data.idx_momp(iW).Edge(j)<=Process_data.idx_momp(iW).track_lost(j)) || ...
                    isnan(Process_data.idx_momp(iW).track_lost(j))
                goodtrack(j) = true;
            else
                ErrorType(j) = 5;
                warning('Unexpected case')
            end
        else            % momp recoreded too early
            assert(Process_data.idx_momp(iW).track_lost(j)<minlength)
            ErrorType(j) = 4;
        end
    end
    assert(sum(goodtrack)==sum(ErrorType==0))


    temp_plot='y';

    for j=1:length(goodtrack)
        if plotting && j<abs(plotting) && temp_plot=='y' && ErrorType(j)
            figure(994);clf
            subplot(121)
            hold on
            plot(T2, FRET(:,j),'c','linewidth',1)
            if Process_data.idx_momp(iW).Edge(j)>0
                plot(Process_data.idx_momp(iW).Edge(j)*[1 1], [0 max(FRET(:))], '.b-')
            end
            title(Process_data.Legends.ErrorType(ErrorType(j)))
            ylims = ylim;

            plot([1 1]*minT, ylims,'-k')
            plot([1 1]*maxT, ylims,'-k')
            xlim([0 max(T2)+10])

            subplot(122)
            hold on
            plot(T2, [0;diff(FRET(:,j))],'k','linewidth',1)
            xlim([0 max(T2)+10])

            if plotting>0
                drawnow
                temp_plot = input('    Do you want more? y/n [y]: ', 's');
                if isempty(temp_plot);temp_plot = 'y';
                elseif length(temp_plot)>1;temp_plot = 'n';end
            else pause(.2);end
        end
    end

    if ~strcmp(Process_data.RFPchannel,'None')
        if isfield(Process_data.rawdata,'rfp0');
            assert(length(Process_data.rawdata(iW).rfp0)==length(goodtrack))
            Process_data.rfp0{iW} = Process_data.rawdata(iW).rfp0(goodtrack);
        elseif isfield(Process_data.rawdata,'rfp') && ...
                isvector(Process_data.rawdata(iW).rfp);
            assert(length(Process_data.rawdata(iW).rfp)==length(goodtrack))
            Process_data.rfp0{iW} = Process_data.rawdata(iW).rfp(goodtrack);
        end
    end
    goodtrack = find(goodtrack);
    Process_data.NgoodTraj(iW) = length(goodtrack);
    Process_data.Traj(iW).goodtrack = goodtrack;
    Process_data.Traj(iW).ErrorType = ErrorType;
    % create empty matrices
    Process_data.goodTraj(iW).FRET = NaN(cutoff2,length(goodtrack));
    Process_data.goodTraj(iW).ICTraj = NaN(cutoff2,length(goodtrack));
    Process_data.goodTraj(iW).C8Activity = NaN(cutoff,length(goodtrack));

    if strcmp(Process_data.RFPchannel,'FLIP')
        assert(all(cell2mat( cellfun2(@(x) any(~isnan(x)), {Process_data.rawdata.rfp}))), ...
            'some trajectories have NaN in the rfp variable whereas FLIP was specified to be used')
        Process_data.goodTraj(iW).FLIPlevel{iW} = max(Process_data.rawdata(iW).rfp(:,goodtrack),1);
        Process_data.goodTraj(iW).FLIPlevel = NaN(cutoff2,length(goodtrack));
        warning('wip: check proper FLIP import')
        ProcessFLIP = true;
    else
        ProcessFLIP = false;
    end
    if isfield(Process_data.rawdata,'cfp')
        Process_data.goodTraj(iW).cfp = NaN(cutoff2,length(goodtrack));
        Process_cfp = true;
    else
        Process_cfp = false;
    end

    Process_data.stats(iW).Surviving = isinf(idx_momp(goodtrack))';
    Process_data.stats(iW).Endidx = NaN(1,length(goodtrack));
    Process_data.stats(iW).MompTime = NaN(1,length(goodtrack));
    Process_data.stats(iW).MompIdx = Process_data.stats(iW).MompTime;
    Process_data.stats(iW).MompThres = Process_data.stats(iW).MompTime;
    Process_data.stats(iW).MaxIC = Process_data.stats(iW).MompTime;
    Process_data.stats(iW).MaxTimeIC = Process_data.stats(iW).MompTime;

    Process_data.stats(iW).MompActivity = Process_data.stats(iW).MompTime;
    Process_data.stats(iW).MaxActivity = Process_data.stats(iW).MompTime;
    Process_data.stats(iW).MaxActivityIdx = Process_data.stats(iW).MompTime;
    Process_data.stats(iW).MaxActivityTime = Process_data.stats(iW).MompTime;

    endidx = NaN(length(goodtrack),1);

    fprintf('\tTrajectories analysis\n')
    temp_plot = 'y';
    for j=1:Process_data.NgoodTraj(iW)
        % loop through the different cells

        Cidx = goodtrack(j);

        % end of the trajectory (last point);
        if isinf(idx_momp(Cidx))  % surviving cell
            if any(isnan(FRET(1:cutoff2,Cidx))) % tracking lost between cutoff and cutoff2
                idx=cutoff;
                nanidx = find(isnan(FRET(:,Cidx)),1,'first');
                maxidx = min(cutoff2, nanidx+ceil(windowSize/3)-1);
                % patch a few representative frames to have a good smoothing
                FRET(nanidx:maxidx,Cidx) = FRET((nanidx-1):-1:(2*nanidx-maxidx-1),Cidx);
                idx2=min(cutoff2, find(~isnan(FRET(:,Cidx)),1,'last'));
            else
                idx=cutoff;
                idx2=cutoff2;
            end
        else        % MOMP
            assert(idx_momp(Cidx)>0)
            idx = min( [find(T2 <= (T2(idx_momp(Cidx)+1)+TafterMOMP),1,'last'), ...
                find(~isnan(FRET(:,Cidx)),1,'last'), cutoff] );
            idx2 = min( [windowSize+find(T2 <= (T2(idx_momp(Cidx)+1)+TafterMOMP),1,'last'), ...
                find(~isnan(FRET(:,Cidx)),1,'last'), cutoff2] );
        end
        endidx(j) = idx;
        Process_data.stats(iW).Endidx(j) = idx;

        % smoothing of the trajectory, remove the control
        if idx2>(3*windowSize)
            Process_data.goodTraj(iW).ICTraj(1:idx2,j) = filtfilt(ones(1,windowSize)/windowSize,1, ...
                FRET(1:idx2,Cidx));
            if ProcessFLIP
                Process_data.goodTraj(iW).FLIPlevel(1:idx2,j) = filtfilt(ones(1,windowSize)/windowSize,1, ...
                    Process_data.rawdata(iW).rfp(1:idx2,Cidx));
            end
            if Process_cfp
                Process_data.goodTraj(iW).cfp(1:idx2,j) = filtfilt(ones(1,windowSize)/windowSize,1, ...
                    Process_data.rawdata(iW).cfp(1:idx2,Cidx));
            end

        else
            Process_data.goodTraj(iW).ICTraj(1:idx2,j) = ...
                smooth(FRET(1:idx2,Cidx), windowSize, 'rlowess');
            if ProcessFLIP
                Process_data.goodTraj(iW).FLIPlevel(1:idx2,j) = ...
                    smooth(Process_data.rawdata(iW).rfp(1:idx2,Cidx), windowSize, 'rlowess');
            end
            if Process_cfp
                Process_data.goodTraj(iW).cfp(1:idx2,j) = ...
                    smooth(Process_data.rawdata(iW).cfp(1:idx2,Cidx), windowSize, 'rlowess');
            end
        end
        % robust minimum: bottom 5%
        robustmin = quantile(Process_data.goodTraj(iW).ICTraj(1:idx,j), .05);
        % normalization by minimum
        Process_data.goodTraj(iW).ICTraj(:,j) = Process_data.goodTraj(iW).ICTraj(:,j) -robustmin;
        Process_data.goodTraj(iW).ICTraj(Process_data.goodTraj(iW).ICTraj(:,j)<1e-5,j) = 1e-5;
        Process_data.goodTraj(iW).FRET(:,j) = FRET(:,Cidx) -robustmin;

        Stdnoise = std(Process_data.goodTraj(iW).ICTraj(1:idx2,j)...
            -Process_data.goodTraj(iW).FRET(1:idx2,j));
        Process_data.goodTraj(iW).FRET(Process_data.goodTraj(iW).ICTraj(:,j)<Stdnoise,j) = ...
            Process_data.goodTraj(iW).ICTraj(Process_data.goodTraj(iW).ICTraj(:,j)<Stdnoise,j);


        if idx_momp(Cidx)>0 &&  idx_momp(Cidx)<cutoff % cell dies before cutoff
            Process_data.stats(iW).MompIdx(j) = idx_momp(Cidx);
            Process_data.stats(iW).MompTime(j) = Process_data.T(idx_momp(Cidx));
            Process_data.stats(iW).MompThres(j) = Process_data.goodTraj(iW).ICTraj(idx_momp(Cidx),j);
        else
            Process_data.stats(iW).MompIdx(j) = Inf;
            Process_data.stats(iW).MompTime(j) = Inf;
        end

        % Evaluate C8 activity (derivate)
        if idx2>(3*windowSize)
            temp = filtfilt(ones(1,windowSize)/windowSize,1, ...
                [0;diff(Process_data.goodTraj(iW).ICTraj(1:min(end,idx2),j))./diff(T2(1:min(end,idx2))')]);
            Process_data.goodTraj(iW).C8Activity(1:idx,j) = temp(1:idx);
        else
            Process_data.goodTraj(iW).C8Activity(1:idx,j) = ...
                smooth([0;diff(Process_data.goodTraj(iW).ICTraj(1:idx2,j))./diff(T2(1:idx2)')],ceil(windowSize/2));
        end


        % removing the end of trajectory
        Process_data.goodTraj(iW).ICTraj((idx+1):end,j) = NaN;
        Process_data.goodTraj(iW).FRET((idx+1):end,j) = NaN;
        if ProcessFLIP
            Process_data.goodTraj(iW).FLIPlevel((idx+1):end,j) = NaN;
        end
        if Process_cfp
            Process_data.goodTraj(iW).cfp((idx+1):end,j) = NaN;
        end

        % record max IC and position
        [Process_data.stats(iW).MaxIC(j), Process_data.stats(iW).MaxTimeIC(j)] = ...
            nanmax(Process_data.goodTraj(iW).ICTraj(:,j));
        Process_data.stats(iW).MaxTimeIC(j) = Process_data.T(Process_data.stats(iW).MaxTimeIC(j));

        % determine max activity
        % temp = max(1e-10,smooth(prop.C8Activity{k,i}(1:min(idx,cutoff-prop.windowSize),j),3));
        % previous code used a smoothing function before determining the min
        % --> replace bin the min removing the last 2 points (to avoid
        % overestimating due to feedback
        temp = max(1e-10,Process_data.goodTraj(iW).C8Activity(1:(idx-2),j));
        [Process_data.stats(iW).MaxActivity(j), Process_data.stats(iW).MaxActivityIdx(j)] = ...
            nanmax(temp);
        Process_data.stats(iW).MaxActivityTime(j) = ...
            Process_data.T(Process_data.stats(iW).MaxActivityIdx(j));


        % plotting
        if plotting && j<abs(plotting) && temp_plot=='y'
            figure(994);clf
            subplot(121)
            hold on
            plot(T2, Process_data.goodTraj(iW).FRET(:,j),'c','linewidth',1)
            plot(T2(1:idx), Process_data.goodTraj(iW).ICTraj(1:idx,j),'m','linewidth',2)
            plot(Process_data.stats(iW).MaxTimeIC(j), Process_data.stats(iW).MaxIC(j), '.k')
            plot(Process_data.stats(iW).MaxActivityTime(j)*[1 1], ...
                [0 Process_data.goodTraj(iW).ICTraj(Process_data.stats(iW).MaxActivityIdx(j),j)], ...
                'ok:','linewidth',2)
            if idx_momp(Cidx)>0
                plot(Process_data.stats(iW).MompTime(j)*[1 1], ...
                    [0 Process_data.stats(iW).MompThres(j)], '.b-')
            end

            subplot(122)
            hold on
            plot(T2(1:idx), [0;diff(Process_data.goodTraj(iW).ICTraj(1:idx,j))./diff(Process_data.T(1:idx)')],...
                'k','linewidth',1)
            plot(T2(1:idx), Process_data.goodTraj(iW).C8Activity(1:idx,j),'b','linewidth',2)
            plot(T2([1 idx]), [0 0],'c')
            plot(Process_data.stats(iW).MaxActivityTime(j)*[1 1], ...
                [0 Process_data.stats(iW).MaxActivity(j)], 'ok-')
            if ~isnan(Process_data.stats(iW).MompIdx(j))
                plot(Process_data.stats(iW).MompTime(j)*[1 1], ...
                    [0 Process_data.stats(iW).MompActivity(j)], '-r.')
            end
            ylim([-.002 .02])
            if plotting>0
                drawnow
                temp_plot = input('    Do you want more? y/n [y]: ', 's');
                if isempty(temp_plot)
                    temp_plot = 'y';
                elseif length(temp_plot)>1
                    temp_plot = 'n';
                end
            else
                pause(.5)
            end
        end

    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% record the MOMP as the max FRET activity (for cells that die)
    %%%%% override the MOMP determined by the edge in further analyses


    Process_data.stats(iW).FRETMompTime = Inf(1,length(goodtrack));
    Process_data.stats(iW).FRETSurviving = true(1,length(goodtrack));
    Process_data.stats(iW).FRETMompTimeIdx = Inf(1,length(goodtrack));
    Process_data.stats(iW).FRETMompActivity = NaN(1,length(goodtrack));
    Process_data.stats(iW).FRETPreMompTimeIdx = Inf(1,length(goodtrack));
    Process_data.stats(iW).FRETPreMompTime = Inf(1,length(goodtrack));
    Process_data.stats(iW).MaxActivityPreFRETMomp = NaN(1,length(goodtrack));
    Process_data.stats(iW).MaxFRETPreFRETMomp = NaN(1,length(goodtrack));
    Process_data.stats(iW).MaxActivityIdxPreFRETMomp = NaN(1,length(goodtrack));

    fprintf('\tFinding the MOMP according to FRET peak activation\n');

    temp_plot = 'y';
    for j = find(~Process_data.stats(iW).Surviving)
        % identify MOMP time based on raw FRET data

        Stdnoise = std(Process_data.goodTraj(iW).ICTraj(1:endidx(j),j)...
            -Process_data.goodTraj(iW).FRET(1:endidx(j),j));

        rawActivity = [0;diff(Process_data.goodTraj(iW).ICTraj(1:endidx(j),j))./...
            diff(Process_data.T(1:endidx(j))')];
        rawDelta = rawActivity - Process_data.goodTraj(iW).C8Activity(1:endidx(j),j);  % delta from smoothed curve
        rawDelta = rawDelta .* ...
            (Process_data.goodTraj(iW).C8Activity(1:endidx(j),j) > ...
            .7*max(Process_data.goodTraj(iW).C8Activity(1:endidx(j),j)) ) ...
            .* (Process_data.goodTraj(iW).ICTraj(1:endidx(j),j)>Stdnoise) .* ...
            ([1:endidx(j)]'>windowSize);
        rawDelta(rawDelta<0) = 0;
        if endidx(j)==cutoff
            trawDelta = rawDelta(1:(end-windowSize));
        else
            trawDelta = rawDelta;
        end

        % identify 5 consecutive frames with increase above 2*std
        HighActivity = (trawDelta>(3*std(trawDelta))) + (trawDelta>(2*std(trawDelta))) - ...
            (trawDelta<(-2*std(trawDelta)));
        HighActivityValues = HighActivity.*sqrt(trawDelta);
        CumulativeHighActivity = MovingWindow(HighActivity,7,@nansum);
        CumulativeHighActivityValues = MovingWindow(HighActivityValues,7,@nansum);
        MompCandidates = find( CumulativeHighActivity>3 );

        if ~isempty(MompCandidates)
            [~,TopMompCandidate] = max( CumulativeHighActivityValues );
            MompCandidates = unique([MompCandidates; TopMompCandidate]);
            MompCandidateValues = CumulativeHighActivity(MompCandidates)/max( CumulativeHighActivity );

            % first identify the burst points (5 frames around Momp)
            BurstIdx = MompCandidates(abs(MompCandidates-TopMompCandidate)<6);
            % find the max activity in the burst (MOMP time):
            [Process_data.stats(iW).FRETMompActivity(j), temp] = max(rawActivity(BurstIdx));
            Process_data.stats(iW).FRETMompTimeIdx(j) = BurstIdx(temp);
            Process_data.stats(iW).FRETMompTime(j) = ...
                Process_data.T(Process_data.stats(iW).FRETMompTimeIdx(j));

            Process_data.stats(iW).FRETPreMompTimeIdx(j) = min(BurstIdx) - SubstractedFrames;
            Process_data.stats(iW).FRETPreMompTime(j) = ...
                Process_data.T( Process_data.stats(iW).FRETPreMompTimeIdx(j) );

            % finally find the last values value prior to the MOMP
            lastframe = min(Process_data.stats(iW).FRETPreMompTimeIdx(j)+2, ...
                find(~isnan(Process_data.goodTraj(iW).C8Activity),1,'last'));
            temp = max(1e-10,smooth(Process_data.goodTraj(iW).C8Activity(1:lastframe,j),3)); % add 2 frames for smoothing
            [Process_data.stats(iW).MaxActivityPreFRETMomp(j), ...
                Process_data.stats(iW).MaxActivityIdxPreFRETMomp(j)] = nanmax(temp(1:(end-2))); % remove 2 frames of the smoothing
            Process_data.stats(iW).MaxFRETPreFRETMomp(j) = Process_data.goodTraj(iW).ICTraj(Process_data.stats(iW).MaxActivityIdxPreFRETMomp(j) ,j);
            assert(~isnan(Process_data.stats(iW).MaxActivityPreFRETMomp(j)))
        else
            MompCandidateValues = [];
        end

        if ~plotting
            assert(all(isnan(Process_data.stats(iW).MaxFRETPreFRETMomp)==...
                isnan(Process_data.stats(iW).FRETMompActivity)))
        end

        if plotting && j<abs(plotting) && temp_plot=='y'
            figure(994);clf
            subplot(121)
            hold on
            plot(T2, Process_data.goodTraj(iW).FRET(:,j),'c','linewidth',1)
            plot(Process_data.T(1:endidx(j)), Process_data.goodTraj(iW).ICTraj(1:endidx(j),j),'m','linewidth',2)
            plot(Process_data.stats(iW).MaxTimeIC(j),Process_data.stats(iW).MaxIC(j), '.k')
            plot(Process_data.stats(iW).MaxActivityTime(j)*[1 1], ...
                [0 Process_data.goodTraj(iW).ICTraj(Process_data.stats(iW).MaxActivityIdx(j),j)], 'ok:','linewidth',2)
            if idx_momp(Cidx)>0
                plot(Process_data.stats(iW).MompTime(j)*[1 1], ...
                    [0 Process_data.stats(iW).MompThres(j)], '.b-')
            end
            for iM=1:length(MompCandidates)
                plot(Process_data.T(MompCandidates(iM))*[1 1],...
                    MompCandidateValues(iM)*[0 max(Process_data.goodTraj(iW).ICTraj(:,j))], 'og-')
            end
            plot(Process_data.stats(iW).FRETMompTime(j)*[1 1], ...
                [0 max(Process_data.goodTraj(iW).ICTraj(:,j))], 'or-')
            plot(Process_data.stats(iW).FRETPreMompTime(j)*[1 1], ...
                [0 max(Process_data.goodTraj(iW).ICTraj(:,j))], 'ok-')

            subplot(122)
            hold on
            plot(Process_data.T(1:endidx(j)), rawActivity(1:endidx(j)),'k','linewidth',1)
            plot(Process_data.T(1:endidx(j)), Process_data.goodTraj(iW).C8Activity(1:endidx(j),j),'b','linewidth',2)
            plot(Process_data.T([1 endidx(j)]), [0 0],'c')
            plot(Process_data.stats(iW).MaxActivityTime(j)*[1 1], ...
                [0 Process_data.stats(iW).MaxActivity(j)], 'ok-')

            for iM=1:length(MompCandidates)
                plot(Process_data.T(MompCandidates(iM))*[1 1], ...
                    MompCandidateValues(iM)*[0 Process_data.stats(iW).MaxActivity(j)], 'og-')
            end
            plot(Process_data.stats(iW).FRETMompTime(j)*[1 1], ...
                [0 Process_data.stats(iW).FRETMompActivity(j)], 'or-')
            plot(Process_data.stats(iW).FRETMompTime(j)*[0 1], ...
                [1 1]*Process_data.stats(iW).FRETMompActivity(j), 'or-')
            plot(Process_data.stats(iW).FRETPreMompTime(j)*[1 1], ...
                [0 Process_data.stats(iW).MaxActivity(j)], 'ok-')
            plot(Process_data.stats(iW).FRETPreMompTime(j)*[0 1], ...
                Process_data.stats(iW).MaxActivityPreFRETMomp(j)*[1 1], 'ok-')


            if ~isnan(Process_data.stats(iW).MompIdx(j))
                plot(Process_data.stats(iW).MompTime(j)*[1 1], ...
                    [0 Process_data.stats(iW).MompActivity(j)], '-r.')
            end
            ylim([-.002 .02])


            if plotting>0
                temp_plot = input('    Do you want more? y/n [y]: ', 's');
                if isempty(temp_plot)
                    temp_plot = 'y';
                elseif length(temp_plot)>1
                    temp_plot = 'n';
                end
            else
                pause(.5)
            end
        end
    end

    % Prior was using: Process_data.stats(iW).FRETMompActivity (note of
    % March 20)
    Process_data.stats(iW).FRETSurviving = isnan(Process_data.stats(iW).MaxActivityPreFRETMomp);
    Process_data.stats(iW).preMompMaxActivity = Process_data.stats(iW).MaxActivityPreFRETMomp;
    Process_data.stats(iW).preMompMaxActivity(Process_data.stats(iW).FRETSurviving) = ...
        Process_data.stats(iW).MaxActivity(Process_data.stats(iW).FRETSurviving);

    Process_data.stats(iW).MaxFRETPreFRETMomp = Process_data.stats(iW).MaxFRETPreFRETMomp;
    Process_data.stats(iW).MaxFRETPreFRETMomp(Process_data.stats(iW).FRETSurviving) = ...
        Process_data.stats(iW).MaxIC(Process_data.stats(iW).FRETSurviving);

    Process_data.stats(iW).preMompMaxActivityIdx = Process_data.stats(iW).MaxActivityIdxPreFRETMomp;
    Process_data.stats(iW).preMompMaxActivityIdx(Process_data.stats(iW).FRETSurviving) = ...
        Process_data.stats(iW).MaxActivityIdx(Process_data.stats(iW).FRETSurviving);

    Process_data.stats(iW).preMompMaxActivityTime = ...
        Process_data.T(Process_data.stats(iW).preMompMaxActivityIdx);


    % removing the end of trajectory
    Process_data.goodTraj(iW).ICTraj((cutoff+1):end,:) = [];
    Process_data.goodTraj(iW).FRET((cutoff+1):end,:) = [];
    if ProcessFLIP
        Process_data.goodTraj(iW).FLIPlevel((cutoff+1):end,:) = [];
    end
    if Process_cfp
        Process_data.goodTraj(iW).cfp((cutoff+1):end,:) = [];
    end

    fprintf(['%i good trajectories (%i survive until %.0fh; %i MOMP) out of %i:\n' ...
        '\t%i bad tracking, %i not tracked until 5h (no MOMP), %i not tracked until %.0fh (no MOMP),\n' ...
        '\t%i unclassified, %i MOMP/track lost before %.1f min\n'],...
        length(Process_data.stats(iW).Surviving), sum(Process_data.stats(iW).Surviving), ...
        maxT/60, sum(~Process_data.stats(iW).Surviving), length(ErrorType), ...
        sum(ErrorType==1), sum(ErrorType==2), sum(ErrorType==12), maxT/60, ...
        sum(ErrorType==3), sum(ErrorType==4), Process_data.T(minlength));


end
