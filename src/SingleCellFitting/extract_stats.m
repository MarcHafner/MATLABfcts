function [stats, allActThreshold, allf_mod, uniqueDoses, Dates] = extract_stats(Exp_data, ActThreshold)
% [stats, allActThreshold, allf_mod, Doses, Date] = extract_stats(Exp_data)
Generate_TRIALplotProp

Drugidx = ismember(Exp_data.ExpKey.Properties.VariableNames, {'TRAIL' 'Mapa' 'Apomab' 'Fc'});
Doses =table2array(Exp_data.ExpKey(:,Drugidx));
uniqueDoses = unique(Doses,'rows');
Dates = unique(Exp_data.ExpKey.Date);

regularterm = 10;

stats = struct;
allSurval = [];
allMOMPval = [];
for iD = 1:size(uniqueDoses,1)

    cnt = 0;
    for iDt = 1:length(Dates)
        idx = find(all(Doses==(ones(size(Doses,1),1)*uniqueDoses(iD,:)),2) &...
            strcmp(Exp_data.ExpKey.Date, Dates(iDt)));

        for i=1:length(idx)
            cnt=cnt+1;

            temp = Exp_data.fits(idx(i)).k;
            temp(temp==0) = 1e-7;
            stats(iD,cnt).k = temp;
            stats(iD,cnt).medk = 10.^(nanmean(log10(temp)));
            stats(iD,cnt).tau = Exp_data.Fate(idx(i)).preMompMaxActivityTime;
            stats(iD,cnt).Surv = mean(Exp_data.Fate(idx(i)).Surviving);


            Surval  = Exp_data.Fate(idx(i)).preMompMaxActivity(:,Exp_data.Fate(idx(i)).Surviving);
            MOMPval = Exp_data.Fate(idx(i)).preMompMaxActivity(:,~Exp_data.Fate(idx(i)).Surviving);
            allSurval = [allSurval Surval];
            allMOMPval = [allMOMPval MOMPval];

            optimThreshold = @(x) sum(x-MOMPval(MOMPval<x)) + sum(Surval(Surval>x)-x);
%             optimThreshold = @(x) mean([MOMPval<x Surval>x]) + ...
%                 regularterm*(sum(x-MOMPval(MOMPval<x)) + sum(Surval(Surval>x)-x));
            stats(iD,cnt).ActThreshold = fminsearch(optimThreshold, 2e-3);
            stats(iD,cnt).Acc = mean([MOMPval>stats(iD,cnt).ActThreshold ...
                Surval<stats(iD,cnt).ActThreshold]);


            tM = Exp_data.Fate(idx(i)).preMompMaxActivityTime;
            stats(iD,cnt).f_mod = @(x) 20+(.5*stats(iD,cnt).ActThreshold)./x;
            stats(iD,cnt).mod_Tacc = mean((tM<stats(iD,cnt).f_mod(temp)) == Exp_data.Fate(idx(i)).Surviving);

        end
    end
end

allActThreshold = [];
allf_mod = [];

if ~exist('ActThreshold','var') || isempty(ActThreshold)
        optimThreshold = @(x) sum(x-allMOMPval(allMOMPval<x)) + sum(allSurval(allSurval>x)-x);
%     optimThreshold = @(x) mean([allMOMPval<x allSurval>x]) + ...
%         regularterm*(sum(x-allMOMPval(allMOMPval<x)) + sum(allSurval(allSurval>x)-x));
    allActThreshold = fminsearch(optimThreshold, 2e-3);
    ActThreshold = allActThreshold;
    allf_mod = @(x) 20+(.5*allActThreshold)./x;
    f_mod = allf_mod;

elseif ActThreshold==-1
    ActThreshold = TRAIL_prop.allActThreshold;
    f_mod = TRAIL_prop.allf_mod;
    allActThreshold = ActThreshold;
else
    f_mod = @(x) 20+(.5*ActThreshold)./x;
end

for iD = 1:size(uniqueDoses,1)

    cnt = 0;
    for iDt = 1:length(Dates)
        idx = find(all(Doses==(ones(size(Doses,1),1)*uniqueDoses(iD,:)),2) &...
            strcmp(Exp_data.ExpKey.Date, Dates(iDt)));

        for i=1:length(idx)
            cnt=cnt+1;
            Surval  = Exp_data.Fate(idx(i)).preMompMaxActivity(:,Exp_data.Fate(idx(i)).Surviving);
            MOMPval = Exp_data.Fate(idx(i)).preMompMaxActivity(:,~Exp_data.Fate(idx(i)).Surviving);
            tM = Exp_data.Fate(idx(i)).preMompMaxActivityTime;
            temp = Exp_data.fits(idx(i)).k;
            temp(temp==0) = 1e-7;

            stats(iD,cnt).allAcc = mean([MOMPval>ActThreshold ...
                Surval<ActThreshold]);
            stats(iD,cnt).allmod_Tacc = mean((tM<f_mod(temp)) == Exp_data.Fate(idx(i)).Surviving);
        end
    end
end
