
function [allGIs, t_rfits] = AngleSynergyAnalysis(t_data, extraGIs, plotting, verbatim)
% [allGIs, t_rfits] = AngleSynergyAnalysis(t_data, extraGIs, plotting, verbatim)
%
% Inputs:
%   t_data:     table with the columns: DrugName, Conc, DrugName2, Conc2,
%               RelGrowth and optional column CellLine (for figure name)
%   extraGIs:   additional GI levels. Default is [25 65 80 90].
%   plotting:   plot nothing (0), only the final results (1, default) or
%                   the fits and the final result (2)
%   verbatim:   disply fit informations (default: true)
%
% Outputs:
%   allGIs:     all results for each GI levels
%   t_rfits:    fits along the ratio of drug concentrations
%
%
% Improvment to make:
%   - support RelCount as an alternative to RelGrowth
%   - GIx is set to 1.5-fold higher than max dose if GIx cannot be defined
%   - normalization currently by GI50
%   - support other input formats
%   - currently assuming that all doses used in the combo are also in the
%       single agent
%
%


if ~exist('plotting','var') || plotting
    global Plotting_parameters
    if ~exist('plotting','var') || islogical(plotting) || plotting<2
        plotting = 1;
    end
else
    plotting = 0;
end

if ~exist('verbatim','var')
    verbatim = true;
end
%%%%%%%%%%%%%%
%%%% assume outputs from the D300 DrugResponse Suite (like t_mean)
%%%%%%%%%%%%%%
if ~exist('extraGIs','var') || isempty(extraGIs)
    extraGIs = [20 35 65 80 95];
else
    extraGIs = unique(ToRow(extraGIs),'stable');
end
if ~ismember(50, extraGIs) % put the GI50 in front
    GIvalues = [50 extraGIs];
    GI50idx = 1;
else
    GIvalues = extraGIs; % keep the order
    extraGIs = setdiff(extraGIs,50,'stable');
    GI50idx = find(GIvalues==50);
end


allGIs = struct('GIval', num2cell(GIvalues), 'BAratio', [], 'RealConc', [], ...
    'NullConc', [], 'CI', [], 'RelBAratio', [], 'CappedGI', [], ...
    'GIconc', [], 'meanCI', [], 'intrpBAratio', [], 'intrpConc', [], ...
    'intrpCI', [], 'intrpRelBAratio', []);

if ~exist('t_data','var') || isempty(t_data)
    if verbatim
        warnprintf('Empty data; returning empty structure'), end
    return
end


if isvariable(t_data,'CellLine')
    assert(length(unique(t_data.CellLine))==1, 'There are more than one cell line: %s', ...
        strjoin(cellstr(unique(t_data.CellLine)), ', '))
end

Drugs = unique(t_data.DrugName(t_data.DrugName2~='-'));
assert(~isempty(Drugs), 'No ''Drug1'' found')
assert(length(Drugs)==1, 'There is more than one ''Drug1'': %s', ...
    strjoin(cellstr(Drugs),', '))

temp = setdiff(unique(t_data.DrugName(t_data.DrugName2=='-')), Drugs);
assert(length(temp)==1, 'There is more than one ''Drug2'': %s', ...
    strjoin(cellstr(temp),', '))
Drugs(2) = temp;

%%%%%%%%%%%%%%
%%%% currently assume growth data with 'RelGrowth'
%%%% do NOT handle other metics
%%%%%%%%%%%%%%

%% evaluate the GI50 for each drug independently

GIs = Inf(length(GIvalues),2); % already set to max concentration
cappedGI = false(length(GIvalues),2);
Hills = NaN(1,2);
Einfs = Hills;
EC50s = Hills;
Emaxs = Hills;
r2s = Hills;
NormFactor = Hills;

if plotting==2
    get_newfigure(999,[460 100 800 700]);
    get_subaxes(2,2,1,[],1)
    set(gca,'xscale','log')

    colors = [.7 .1 0; .2 0 1];
end

for i=1:2
    t_sub = sortrows(t_data(t_data.DrugName==Drugs(i) & t_data.DrugName2=='-',:), ...
        'Conc');
    Emaxs(i) = 1-min(t_sub.RelGrowth);

    opt = struct('pcutoff',.5,'extrapolrange',10^1.5,'Robust',true);
    [GIs(GI50idx,i), Hills(i), Einfs(i), ~, r2s(i), EC50s(i), fitfct] = ...
        ICcurve_fit(t_sub.Conc, t_sub.RelGrowth, 'GI50', opt);
    %%%% replace by the actual subfunction and stack it in the file
    %%%%%%%%%%
    if verbatim
        fprintf('\tfit for %-12s: EC50=%.3g, GI50=%.3g, r2=%.2f\n', char(Drugs(i)), ...
            EC50s(i), GIs(GI50idx,i), r2s(i)); end
    if verbatim
        warnassert(r2s(i)>.7, 'Bad fit for %s !', char(Drugs(i))), end
    %%%%% should that stop the algorithm ??
    if ~isfinite(EC50s(i))
        if verbatim
            warnprintf('data for drug %s cannot be fitted --> all GIsets to min/max dose', ...
                char(Drugs(i))), end
        cappedGI(:,i) = true;
        GIs((1-(GIvalues/100))<fitfct.b,i) = max(t_sub.Conc)*10^1.5;
        GIs((1-(GIvalues/100))>fitfct.b,i) = min(t_sub.Conc)/10^1.5;
        continue
    end

    % evaluation of the other GIs
    %%%%% repeated code; could be paste in a subfunction
    fitopts = optimoptions('fsolve', 'Display','none');
    for j=1:length(extraGIs)
        if Einfs(i)>extraGIs(j)/100
            GIs(extraGIs(j)==GIvalues,i) = fsolve(@(x) fitfct(x)-(1-extraGIs(j)/100), ...
                EC50s(i), fitopts);
        else
            GIs(extraGIs(j)==GIvalues,i) = Inf;
        end
    end
    %%%%%% the other GI should directly be evaluated in the ICcurve_fit
    %%%%%% function and capped as the GI50 if the values are extrapolated

    cappedGI(:,i) = GIs(:,i)==Inf;
    if any(cappedGI(:,i)) && verbatim
        warnprintf('Values %s are set to maximum for drug %s', ...
            strjoin(num2cellstr(GIvalues(cappedGI(:,i))),','), char(Drugs(i)))
    end
    GIs(cappedGI(:,i),i) = max(t_sub.Conc)*10^1.5;

    if plotting==2
        plot(t_sub.Conc, t_sub.RelGrowth, 'o', 'color', colors(i,:))
        x = 10.^(log10(min(t_sub.Conc)):.01:log10(max(t_sub.Conc)));
        plot(x, fitfct(x), '-', 'color', colors(i,:))
        plot(x(end)*[1 1.2], [1 1]-Emaxs(i), ':', 'color', colors(i,:))
        plot(x(end)*[1 1.2], [1 1]-Einfs(i), '--', 'color', colors(i,:))
        for j=1:length(GIvalues)
            plot(GIs(j,i)*[1 1 1], [0 1-GIvalues(j)/100 1], '.-', 'color', colors(i,:))
        end
    end
end


    %%%%%%%%%%%%%%
    %%%%  how to deal with normalization when there is no GI50 for normalization ????
    %%%%  use EC50??? instead (with a warning?)
    %%%%
    %%%%  currently forcing a GI50 1.5-fold higher than the highest dose
    %%%%  tested !!
    %%%%%%%%%%%%%%

    if all(~cappedGI(GI50idx,:))
        NormConc = 'GI50';
        NormFactor = GIs(GI50idx,:);
    else
        NormGI = max(GIvalues(all(~cappedGI,2)));

        if isempty(NormGI)
            NormConc = 'capGI50';
            NormFactor = GIs(GI50idx,:);
        else
            NormConc = ['GI' num2str(NormGI)];
            NormFactor = GIs(NormGI==GIvalues,:);
            if verbatim
                warnprintf('Normalization by GI%i instead of GI50', NormGI), end
        end
    end

%% caculate the ratios on which to fit the sigmoidal cruves

t_combodata = t_data(t_data.DrugName==Drugs(1) & t_data.DrugName2==Drugs(2),...
    {'Conc' 'Conc2' 'RelGrowth'});

BAratio = 1e-2*round(log10(t_combodata.Conc2./t_combodata.Conc)*1e2);
t_combodata = [t_combodata table(BAratio) ...
    table(t_combodata.Conc+t_combodata.Conc2, 'variablename', {'ConcTot'}) ];
BAratio = unique(t_combodata.BAratio,'sorted');
ratiocounts = hist(t_combodata.BAratio, BAratio);

BAratio = BAratio(ratiocounts>=4);
if verbatim
    fprintf('\tFitting for %i ratios:\n', length(BAratio));end

t_rfits = [table(BAratio) array2table(NaN(length(BAratio),4), ...
    'variablenames', {'Hill' 'Einf' 'r2' 'Emax'})];
colors = parula(length(BAratio));

ratioGIs = Inf(length(BAratio),length(GIvalues));

if plotting==2
    get_subaxes(2,2,2,[],1)
    set(gca,'xscale','log')
end
for i = 1:height(t_rfits)
    t_sub = sortrows(t_combodata(t_combodata.BAratio == t_rfits.BAratio(i),:), ...
        'ConcTot');
    t_rfits.Emax(i) = 1-min(t_sub.RelGrowth);

    %%%%%%%
    opt = struct('pcutoff',1,'extrapolrange',10^1.5);
    [ratioGIs(i,GI50idx), t_rfits.Hill(i), t_rfits.Einf(i), ~, t_rfits.r2(i), ...
        ~, fitfct, ~, log] = ...
        ICcurve_fit(t_sub.ConcTot, t_sub.RelGrowth, 'GI50', opt);
    %%%% replace by the actual subfunction and stack it in the file
    %%%%%%%%%%

    %%%%% repeated code; could be paste in a subfunction
    for j=1:length(extraGIs)
        jidx = extraGIs(j)==GIvalues;
        fitopts = optimoptions('fsolve', 'Display','none');
        temp = fsolve(@(x) fitfct(x)-(1-extraGIs(j)/100), ...
            median(t_sub.ConcTot), fitopts);
        % capping the values to avoid too much extrapolation
        if temp<min(t_sub.ConcTot)*10^-1.5 || imag(temp)~=0
            ratioGIs(i,jidx) = -Inf;
        elseif temp>10^1.5*max(t_sub.ConcTot)
            ratioGIs(i,jidx) = Inf;
        else
            ratioGIs(i,jidx) = temp;
        end
    end
    %%%%%% the other GI should directly be evaluated in the ICcurve_fit
    %%%%%% function and capped as the GI50 if the values are extrapolated
    if verbatim
        fprintf('\t\tfit for BAratio %-5g: GI50=%.3g, r2=%.2f\n', t_rfits.BAratio(i), ...
            ratioGIs(i,GI50idx), t_rfits.r2(i)); end


    if isinf(ratioGIs(i,GI50idx)) || imag(ratioGIs(i,GI50idx))~=0
        if imag(ratioGIs(i,GI50idx))~=0
            ratioGIs(i,GI50idx) = -Inf;
        end
        if verbatim, fprintf(['\t\t --> ' log '\n']), end
    end

    if plotting==2
        plot(t_sub.ConcTot, t_sub.RelGrowth, 'o', 'color', colors(i,:))
        x = 10.^(log10(infmin([ratioGIs(i,GI50idx);t_sub.ConcTot])):.01:...
            log10(infmax([ratioGIs(i,GI50idx);t_sub.ConcTot])));
        plot(x, fitfct(x), '-', 'color', colors(i,:))
        plot(x(end)*[1 1.2], [1 1]-t_rfits.Emax(i), ':', 'color', colors(i,:))
        plot(x(end)*[1 1.2], [1 1]-t_rfits.Einf(i), '--', 'color', colors(i,:))

        plot(ratioGIs(i,:), 1-GIvalues/100, '.', 'color', colors(i,:), ...
            'markersize', 15)
    end

end
t_rfits.RelBAratio = t_rfits.BAratio +log10(NormFactor(1)) -log10(NormFactor(2));

%% calculate the GIs along the directions of sampling (no Hill or Emax here)
Conc = cell(2,1);
AlignConc = cell(2,1);
for i=1:2
    % assuming that all doses used in the combo are also in the single
    % agent
    Conc{i} = unique(t_data.Conc(t_data.DrugName==Drugs(i)),'sorted');
end

alignedGIs = cell(2,1);

for iC=1:2
    if iC==1 % not the most elegant solution
        n = hist_cat(t_combodata.Conc,Conc{iC});
    else
        n = hist_cat(t_combodata.Conc2,Conc{iC});
    end
    idx = find(n>=4);
    AlignConc{iC} = Conc{iC}(idx);
    alignedGIs{iC} = Inf(length(idx), length(GIvalues));

    if verbatim, fprintf('\tFitting for Conc%i (%i):\n', iC, length(idx)), end
    if plotting==2
        colors = parula(length(idx));
        get_subaxes(2,2,2+iC,[],1)
        set(gca,'xscale','log')
    end
    for i=1:length(idx)

        if iC==1 % swapping to keep code the same
            t_sub = t_combodata(t_combodata.Conc==Conc{iC}(idx(i)),:);
        else
            t_sub = t_combodata(t_combodata.Conc2==Conc{iC}(idx(i)),:);
            t_sub.Properties.VariableNames({'Conc' 'Conc2'}) = {'Conc2' 'Conc'};
        end

        RelGrowth0 = mean(t_data.RelGrowth(t_data.DrugName==Drugs(iC) & t_data.DrugName2=='-' & ...
            t_data.Conc==Conc{iC}(idx(i))));

        ranges = [
            RelGrowth0-.05 max([RelGrowth0;t_sub.RelGrowth])+.05  %E0
            0 1    %Emax
            max(min(t_sub.Conc2)*1e-4,1e-7) min(max(t_sub.Conc2)*1e2, 1e3)  %E50
            .1 5    % HS
            ]';
        opt = struct('pcutoff',.9,'extrapolrange',10^1.5,'ranges',ranges, ...
            'Robust', true);
        [alignedGIs{iC}(i,GI50idx), ~, ~, ~, r2, ~, fitfct, ~, log, flag] = ...
            ICcurve_fit(t_sub.Conc2, t_sub.RelGrowth, 'GI50', opt);

        if r2<.7 || flag==0
            if verbatim, fprintf(['\t\t  Bad fit --> ' log '\n']), end
            alignedGIs{iC}(i,:) = NaN;
            continue
        end

        %%%%% repeated code; could be paste in a subfunction
        for j=1:length(extraGIs)
            jidx = extraGIs(j)==GIvalues;
            fitopts = optimoptions('fsolve', 'Display','none');
            [temp, fval, flag] = fsolve(@(x) fitfct(x)-(1-extraGIs(j)/100), ...
                infmax(alignedGIs{iC}(i,GI50idx), fitfct.c), fitopts);
            if imag(temp)~=0
                    error('fsolve yields complex number')
            end
            if flag>=0 % convergence
                % capping the values to avoid too much extrapolation
                if temp<min(t_sub.Conc2)*10^-1.5
                    alignedGIs{iC}(i,jidx) = -Inf;
                elseif temp>10^1.5*max(t_sub.Conc2)
                    alignedGIs{iC}(i,jidx) = Inf;
                else
                    alignedGIs{iC}(i,jidx) = temp;
                end
            elseif min(t_sub.RelGrowth)>(1-extraGIs(j)/100)
                    alignedGIs{iC}(i,jidx) = Inf;
            elseif max(t_sub.RelGrowth)<(1-extraGIs(j)/100)
                    alignedGIs{iC}(i,jidx) = -Inf;
            elseif flag==-2
                if fitfct.b>(1-extraGIs(j)/100) && ...
                        min(t_sub.RelGrowth)<(1-extraGIs(j)/100)
                    alignedGIs{iC}(i,jidx) = ...
                        t_sub.Conc2(argmin(abs(t_sub.RelGrowth-(1-extraGIs(j)/100))));
                elseif fitfct.a<(1-extraGIs(j)/100) && ...
                        max(t_sub.RelGrowth)>(1-extraGIs(j)/100)
                    alignedGIs{iC}(i,jidx) = ...
                        t_sub.Conc2(argmin(abs(t_sub.RelGrowth-(1-extraGIs(j)/100))));
                else
                    error('unknown case for fsolve')
                end
            else
                error('unknown case for fsolve')
            end


        end

        if verbatim
            fprintf('\t\t%-5g: GI50=%.3g, r2=%.2f\n', Conc{iC}(idx(i)), ...
                alignedGIs{iC}(i,GI50idx), r2); end
        if isinf(alignedGIs{iC}(i,GI50idx)) && min(t_sub.RelGrowth)>.55 && ...
                max(t_sub.RelGrowth)<.45
            if verbatim, fprintf(['\t\t  --> ' log '\n']), end
        end
        if imag(alignedGIs{iC}(i,GI50idx))~=0
            error('issue with fitting')
        end


        if plotting==2
            plot(t_sub.Conc2, t_sub.RelGrowth, 'o', 'color', colors(i,:))

            x = 10.^(log10(infmin([alignedGIs{iC}(i,GI50idx);t_sub.Conc2])):.01:...
                log10(infmax([alignedGIs{iC}(i,GI50idx);t_sub.Conc2])));

            plot(x, fitfct(x), '-', 'color', colors(i,:))
            plot(x(1)*[1 1.2], [1 1]*mean(ranges(:,1)), '--', 'color', colors(i,:))

            plot(alignedGIs{iC}(i,:), 1-GIvalues/100, '.', 'color', colors(i,:), ...
                'markersize', 15)
        end

    end

end



%% getting the actual Concentration for the ratioGIxx
% and stacking with the alignedGIs


w = [ones(height(t_rfits),1) 10.^t_rfits.BAratio];
w = w ./ (sum(w,2)*[1 1]);


% smoothing parameter
Nextra = 4;
Nsmooth = 5;

if plotting==2
    get_newfigure(998,[750 50 900 700])
end
for i=1:length(GIvalues)
    %%
    allGIs(i).CappedGI = cappedGI(i,:);
    allGIs(i).GIconc = GIs(i,:);

    BAratio = real([ t_rfits.BAratio;
        log10(alignedGIs{1}(:,i)./AlignConc{1});
        log10(AlignConc{2}./alignedGIs{2}(:,i))]);
    RealConc = [repmat(ratioGIs(:,i),1,2).*w;
        [AlignConc{1} alignedGIs{1}(:,i)]
        [alignedGIs{2}(:,i) AlignConc{2}]];

    idx = find(isfinite(BAratio) & ...
        BAratio<log10(2*(max(t_combodata.Conc2)/min(t_combodata.Conc))) & ...
        BAratio>log10(2*(min(t_combodata.Conc2)/max(t_combodata.Conc))) & ...
        all(isfinite(RealConc),2));

    NullBAratio = -6:.1:6;
    w2 = [ones(length(NullBAratio),1) 10.^NullBAratio'];
    w2 = w2 ./ (sum(w2,2)*[1 1]);
    allGIs(i).NullConc = w2.*repmat(GIs(i,:),size(w2,1),1);


    if isempty(idx)
        % case where the GIx is not defined.
        if any(isfinite(GIs(i,:)))
            temp = [alignedGIs{1}(:,i); alignedGIs{2}(:,i)];
            if all(isnan(temp))
                DefaultSynergyVal = NaN;
                if verbatim
                    warnprintf('All fits failed for GI%i -> discarded',GIvalues(i)), end
            end
            if all(temp(~isnan(temp))==-Inf)
                % cases of synergy: the GIxx is lower than the combo tested
                DefaultSynergyVal = 1/7;
                if verbatim
                    warnprintf('No GI%i lower than doses in the combo -> synergy',GIvalues(i)), end
            elseif all(temp(~isnan(temp))==Inf)
                DefaultSynergyVal = 7;
                if verbatim
                    warnprintf('GI%i higher than doses in the combo -> antagonist',GIvalues(i)), end
            else
                DefaultSynergyVal = NaN;
                if verbatim
                    warnprintf('No fit within boundaries for GI%i -> discarded',GIvalues(i)), end
            end

            allGIs(i).BAratio = [-1:.1:1];
            allGIs(i).RealConc = NaN(length(allGIs(i).BAratio),2);

            allGIs(i).intrpBAratio = allGIs(i).BAratio;
            allGIs(i).intrpConc = allGIs(i).RealConc;

            allGIs(i).CI = DefaultSynergyVal*ones(length(allGIs(i).BAratio),1);
            % that will plot a dash line
            allGIs(i).intrpCI = DefaultSynergyVal*...
                [1;repmat([1;NaN;1],ceil(length(allGIs(i).BAratio)/3),1);1];
            allGIs(i).intrpCI = allGIs(i).intrpCI(1:length(allGIs(i).BAratio));

            allGIs(i).RelBAratio = allGIs(i).BAratio ...
                +log10(NormFactor(1)) -log10(NormFactor(2));
            allGIs(i).intrpRelBAratio = allGIs(i).RelBAratio;
            allGIs(i).meanCI = DefaultSynergyVal;
        end
        continue
    end


    idx = idx(sortidx(BAratio(idx)));
    allGIs(i).BAratio = BAratio(idx)';

    allGIs(i).RealConc = RealConc(idx,:);
    w2 = [ones(length(allGIs(i).BAratio),1) 10.^allGIs(i).BAratio'];
    w2 = w2 ./ (sum(w2,2)*[1 1]);

    % allGIs(i).CI = sum(allGIs(i).RealConc,2) ./ (w2*GIs(i,:)');
    allGIs(i).CI = sum((allGIs(i).RealConc ./ repmat(GIs(i,:), size(allGIs(i).RealConc,1),1)),2) ./ ...
        sum((w2.*repmat(GIs(i,:),size(w2,1),1))./repmat(GIs(i,:), size(allGIs(i).RealConc,1),1),2);




    % smoothing
    Nextra = 4;
    BAratio = [.1*floor(min(allGIs(i).BAratio)*10)-((Nextra-1):-1:0)*.1-.5 ...
        allGIs(i).BAratio ...
        .1*ceil(max(allGIs(i).BAratio)*10)+(0:(Nextra-1))*.1+.5 ];
    % add a bit of noise to avoid equal ratios (issue for interpolation)
    BAratio = BAratio+1e-5+(rand(size(BAratio))-.5);
    allGIs(i).intrpBAratio = ...
        ((infmin(allGIs(i).BAratio)-.25):.2:(infmax(allGIs(i).BAratio)+.25))';


    c1 = [repmat(GIs(i,1),Nextra,1);  ...
        allGIs(i).RealConc(:,1); ...
        -Inf(Nextra ,1)];
    c1(c1==-Inf) = infmin([t_combodata.Conc;allGIs(i).RealConc(:,1)])*10^-1.5;


    c2 = [-Inf(Nextra ,1); ...
        allGIs(i).RealConc(:,2); ...
        repmat(GIs(i,2),Nextra ,1)];
    c2(c2==-Inf) = infmin([t_combodata.Conc2;allGIs(i).RealConc(:,2)])*10^-1.5;

%     [tempc1,idx1] = unique(c1,'sorted');
%     [tempc2,idx2] = unique(c2,'sorted');
    allGIs(i).intrpConc = 10.^[smooth(interp1(BAratio', log10(c1), ...
        allGIs(i).intrpBAratio, 'PCHIP'),Nsmooth) ...
        smooth(interp1(BAratio, log10(c2), ...
        allGIs(i).intrpBAratio, 'PCHIP'),Nsmooth)];

    CI = [ones(Nextra ,1); allGIs(i).CI; ones(Nextra ,1)];
    allGIs(i).intrpCI = 2.^smooth(interp1(BAratio, log2(CI), ...
        allGIs(i).intrpBAratio, 'PCHIP'),Nsmooth);

    if plotting==2
        % plot QC
        get_subaxes(length(GIvalues),3,i,1,1)
        plot(allGIs(i).BAratio, log10(allGIs(i).RealConc(:,1)), 'x')
        plot(BAratio, log10(c1), '.', 'markersize',12)
        plot(allGIs(i).intrpBAratio, log10(allGIs(i).intrpConc(:,1)), '-')
        plot(NullBAratio, log10(allGIs(i).NullConc(:,1)), '--k')
        plot(-log10(NormFactor(1))+log10(NormFactor(2))*[1 1], [-2 1], '-k')
        plot([0 0], [-2 1], '-', 'color', [.5 .5 .5])
        ylim([-3.5 1.5])

        get_subaxes(length(GIvalues),3,i,2,1)
        plot(allGIs(i).BAratio, log10(allGIs(i).RealConc(:,2)), 'x')
        plot(BAratio, log10(c2), '.', 'markersize',12)
        plot(allGIs(i).intrpBAratio, log10(allGIs(i).intrpConc(:,2)), '-')
        plot(NullBAratio, log10(allGIs(i).NullConc(:,2)), '--k')
        plot(-log10(NormFactor(1))+log10(NormFactor(2))*[1 1], [-2 1], '-k')
        plot([0 0], [-2 1], '-', 'color', [.5 .5 .5])
        ylim([-3.5 1.5])


        get_subaxes(length(GIvalues),3,i,3,1)
        plot(allGIs(i).BAratio, log2(allGIs(i).CI), 'x')
        plot(BAratio, log2(CI), '.', 'markersize',12)
        plot(allGIs(i).intrpBAratio, log2(allGIs(i).intrpCI), '-')
        plot(-log10(NormFactor(1))+log10(NormFactor(2))*[1 1], [-2 1], '-k')
        plot([0 0], [-2 1], '-', 'color', [.5 .5 .5])
        ylim([-3.5 1.5])
    end

    allGIs(i).RelBAratio = allGIs(i).BAratio +log10(NormFactor(1)) -log10(NormFactor(2));
    allGIs(i).intrpRelBAratio = allGIs(i).intrpBAratio +log10(NormFactor(1)) -log10(NormFactor(2));

    % geometric mean of the CI values
    if sum(abs(allGIs(i).intrpRelBAratio)<1)>3
        allGIs(i).meanCI = 2^mean(log2(allGIs(i).intrpCI(abs(allGIs(i).intrpRelBAratio)<1)));
    else
        idx = abs(allGIs(i).intrpRelBAratio - median(allGIs(i).intrpRelBAratio))<1;
        allGIs(i).meanCI = 2^mean(log2(allGIs(i).intrpCI(idx)));
        if verbatim
            warnprintf('mean CI evaluated on RelBA = [%.1f %.1f] for GI%i', ...
                min(allGIs(i).intrpRelBAratio(idx)), max(allGIs(i).intrpRelBAratio(idx)), ...
                GIvalues(i)), end
    end


end

%% smoothing of Hill coeff and Emax/inf
intrpxlims = [.1*floor(min([allGIs.RelBAratio])*10-1) ...
    .1*ceil(max([allGIs.RelBAratio])*10+1)];

RelBAratio = [intrpxlims(1)-(1:(Nextra-1))*.2 ...
    t_rfits.RelBAratio' intrpxlims(2)+(1:(Nextra-1))*.2];
intrpRelrat = (intrpxlims(1)-.25):.2:(intrpxlims(2)+.25);
intrpHill = smooth(interp1(RelBAratio, ...
    [repmat(Hills(1),1,Nextra) t_rfits.Hill' repmat(Hills(2),1,Nextra)], ...
        intrpRelrat, 'PCHIP'),Nsmooth);
intrpEinf = smooth(interp1(RelBAratio, ...
    [repmat(Einfs(1),1,Nextra) t_rfits.Einf' repmat(Einfs(2),1,Nextra)], ...
        intrpRelrat, 'PCHIP'),Nsmooth);
intrpEmax = smooth(interp1(RelBAratio, ...
    [repmat(Emaxs(1),1,Nextra) t_rfits.Emax' repmat(Emaxs(2),1,Nextra)], ...
        intrpRelrat, 'PCHIP'),Nsmooth);



%% plotting the final figure
if plotting==0
    return
end

if isvariable(t_data,'CellLine')
get_newfigure(100,[50 50 700 550], [char(unique(t_data.CellLine)) ...
    '_Combo_' char(Drugs(1)) '+' char(Drugs(2)) '.pdf'])
else
get_newfigure(100,[50 50 700 550], ...
    [ 'Combo_' char(Drugs(1)) '+' char(Drugs(2)) '.pdf'])
end


% spacing for checkboard
xwidth = .02;
ywidth = .03;
space = .005;

% position of axes
pos = [.06 .57 .5 .33; % CI plot
    .67 .67 .3 .23]; % Hill plot
pos(3,:) = [pos(2,1) pos(2,2)-pos(2,4)-2*space pos(2,3:4)]; % Emax plot
pos(5,:) = [pos(2,1) .06 pos(2,3:4)]; % mean CI
pos(4,:) = [pos(1,1) pos(5,2) pos(1,3)-.1 .34+ywidth]; % checkboard

% position for checkboard
subpos = [
    pos(4,1) pos(4,2)+ywidth+space xwidth pos(4,4)-ywidth;  % Drug 1
    pos(4,1)+xwidth+space pos(4,2) pos(4,3)-xwidth ywidth;  % Drug 2
    pos(4,1)+xwidth+space pos(4,2)+ywidth+space pos(4,3)-xwidth pos(4,4)-ywidth; %combo
    pos(4,1)+pos(4,3)+5*space pos(4,2)+.6*pos(4,4) xwidth/2 pos(4,4)*.4; % colorbar
    pos(4,1)+pos(4,3)+2*space pos(4,2) 0 0]; % legend

% define the ticks for the ratio axis
xticks = [intrpxlims(1)-.1 (intrpxlims(1)+.2):(.1*round(diff(intrpxlims)*1.2)):(intrpxlims(2)-.2) ...
    intrpxlims(2)+.1];
xlims = [intrpxlims(1)-.1 intrpxlims(2)+.1];
Concxticks = xticks -log10(NormFactor(1)) +log10(NormFactor(2));
% ticks in ratio of concentration
Concxticklabels = cell(length(xticks),1);
for i=2:length(xticks)-1
    if Concxticks(i)<0
        Concxticklabels{i} = sprintf('%.1f:1', 10^(-Concxticks(i)));
    else
        Concxticklabels{i} = sprintf('1:%.1f', 10^(Concxticks(i)));
    end
end
Concxticklabels{1} = char(Drugs(1));
Concxticklabels{end} = char(Drugs(2));

% ticks in ratio relative to GI50/EC50
Relxticklabels = cell(length(xticks),1);
for i=2:length(xticks)-1
    if xticks(i)<0
        Relxticklabels{i} = sprintf('%.1f:1', 10^(-xticks(i)));
    else
        Relxticklabels{i} = sprintf('1:%.1f', 10^(xticks(i)));
    end
end
Relxticklabels{1} = char(Drugs(1));
Relxticklabels{end} = char(Drugs(2));
xvals = [xticks(1) t_rfits.RelBAratio' xticks(end)];

GIcolors = min([.2+(sum(extraGIs<50):-1:1)'*[0 .2 .8]/sum(extraGIs<50);
    0 0 0;
    .2+(1:sum(extraGIs>50))'*[.4 .8 0]/sum(extraGIs>50)],1);
GIcolors(sortidx(GIvalues,'ascend'),:) = GIcolors;


%%%%%%%%%%%%%%%
% first plot for the CI
% conc ratio axis
get_newaxes(pos(1,:),1)
set(gca,'xtick',xticks,'xticklabel',Concxticklabels,'xticklabelrotation',90, ...
    'ytick',[],'xaxislocation','top',Plotting_parameters.axes{:})
xlim(xlims)
xlabel('Drug ratios (Conc)',Plotting_parameters.axislabel{:})


get_newaxes(pos(1,:),1)
plot(xvals([1 end]), [0 0], '-k')
plot(repmat([xvals([1 end]) NaN],1,2), reshape([1;1;NaN]*[-1 1],[],1)', ...
    '-', 'color', [.7 .7 .7])
plot([0 0], [-32 32], '-k')
plot(log10(NormFactor(1))-log10(NormFactor(2))*[1 1], [-32 32], '-', 'color', [.5 .5 .5])

h = zeros(length(GIvalues),1);
for j=[ToRow(find(GIvalues~=50)) GI50idx] % plot GI50 last
    h(j) = plot(allGIs(j).RelBAratio, log2(allGIs(j).CI), ...
        '.','color',GIcolors(j,:), 'markersize', 12);
    h(j) = plot(allGIs(j).intrpRelBAratio, log2(allGIs(j).intrpCI), ...
        '-','color',GIcolors(j,:),'linewidth',1);
end

set(gca,'xtick',xticks,'xticklabel',Relxticklabels,'xticklabelrotation',90, ...
    'ytick',-5:5,'yticklabel',[strcat('1/',num2cellstr(2.^(5:-1:1))) ...
    strcat(num2cellstr(2.^(0:5)),'/1')],Plotting_parameters.axes{:}, 'box','on')
ylim(log2([min([.45 cellfun(@min,{allGIs(~cellfun(@isempty, {allGIs.BAratio})).CI})])/1.1 ...
    max([2.2 cellfun(@max,{allGIs(~cellfun(@isempty, {allGIs.BAratio})).CI})])*1.1]))
xlim(xlims)
ylabel('Combination index',Plotting_parameters.axislabel{:})
xlabel(['Drug ratios (' NormConc '/' NormConc ' normalized)'],...
    Plotting_parameters.axislabel{:})
% legend(h,[strcat('GI', num2cellstr(GIvalues)) 'Failed points'],'location','best', ...
%     Plotting_parameters.axislabel{:})


%%%%%%%%%%%%%%%%%%
% plot the average CI as a function of GIx
get_newaxes(pos(5,:),1)
plot([0 200], [0 0], '-k')
plot(repmat([0 200 NaN],1,2), reshape([1;1;NaN]*[-1 1],[],1)', ...
    '-', 'color', [.7 .7 .7])

idx = sortidx(GIvalues);
plot(GIvalues(idx), log2([allGIs(idx).meanCI]), '.-k')

set(gca,'ytick',-5:5,'yticklabel',[strcat('1/',num2cellstr(2.^(5:-1:1))) ...
    strcat(num2cellstr(2.^(0:5)),'/1')],Plotting_parameters.axes{:}, 'box','on')
ylim(log2([min([.45 cellfun(@min,{allGIs(~cellfun(@isempty, {allGIs.BAratio})).CI})])/1.1 ...
    max([2.2 cellfun(@max,{allGIs(~cellfun(@isempty, {allGIs.BAratio})).CI})])*1.1]))
xlim([20 max(GIvalues)+10])
ylabel('Combination index',Plotting_parameters.axislabel{:})
xlabel('GI value',Plotting_parameters.axislabel{:})


%%%%%%%%%%%%%%%%%%%
% plot for the Hill coeff
% get_newaxes(pos(2,:),1)
% set(gca,'xtick',xticks,'xticklabel',Concxticklabels,'xticklabelrotation',90, ...
%     'ytick',[],'xaxislocation','top',Plotting_parameters.axes{:})
% xlim(xlims)
get_newaxes(pos(2,:),1)
plot([0 0], [-32 32], '-k')
plot(log10(NormFactor(1))-log10(NormFactor(2))*[1 1], [-32 32], '-', 'color', [.5 .5 .5])
plot(t_rfits.RelBAratio, t_rfits.Hill, '.k')
plot(intrpRelrat, intrpHill , '-k')

xlim(xlims)
ylim([min([.4 t_rfits.Hill']) max([4 t_rfits.Hill'])])
set(gca,'xtick',xticks,'xticklabel',Concxticklabels,'xticklabelrotation',90, ...
    'ytick',0:5,'xaxislocation','top',Plotting_parameters.axes{:}, 'box','on')
xlabel('Drug ratios (Conc)',Plotting_parameters.axislabel{:})
ylabel('Hill coefficient',Plotting_parameters.axislabel{:})

%%%%%%%%%%%%%%%%%%%%%%%
% plot for the Emax/inf
% get_newaxes(pos(3,:),1)
% set(gca,'xtick',xticks,'xticklabel',Concxticklabels,'xticklabelrotation',90, ...
%     'ytick',[],'xaxislocation','top',Plotting_parameters.axes{:})
% xlim(xlims)
get_newaxes(pos(3,:),1)
plot([0 0], [-32 32], '-k')
plot(log10(NormFactor(1))-log10(NormFactor(2))*[1 1], [-32 32], '-', 'color', [.5 .5 .5])
plot(xvals([1 end]), [0 0], '-', 'color', [.7 .7 .7])
h = plot(xvals, 1-[Einfs(1) t_rfits.Einf' Einfs(2)], 'xk');
h(2) = plot(xvals, 1-[Emaxs(1) t_rfits.Emax' Emaxs(2)], '.k');
h = plot(intrpRelrat, 1-intrpEinf, '--k');
h(2) = plot(intrpRelrat, 1-intrpEmax, '-k');
set(gca,'xtick',xticks,'xticklabel',Relxticklabels,'xticklabelrotation',90,...
    Plotting_parameters.axes{:}, 'box','on')
xlim(xlims)
ylim([min([-.3 1-t_rfits.Emax' 1-Emaxs]) max([.3 1-t_rfits.Emax' 1-Emaxs])])
ylabel('Lowest growth',Plotting_parameters.axislabel{:})
xlabel(['Drug ratios (' NormConc '/' NormConc ' normalized)'],...
    Plotting_parameters.axislabel{:})
legend(h,{'1-E_{inf}' '1-E_{max}'},'location','best',...
    Plotting_parameters.axislabel{:})



% plotting the checkboard for combo

%%%%%%% put the following in its own function

%%%%%% assuming carthesian product on conc; need to be generalized
% needs to be split in the columns for single dose and matrix for combos
% (for cases where there is combo are not carthesian product on conc)
%%%
logC1 = log10(Conc{1});
logC2 = log10(Conc{2});
GrowthMx = NaN(length(logC1), length(logC2));
GrowthMx1 = NaN(length(logC1), 1);
GrowthMx2 = NaN(1, length(logC2));

for iC1 = 1:length(logC1)
    GrowthMx1(iC1) = t_data.RelGrowth(t_data.DrugName==Drugs(1) & ...
        t_data.Conc==Conc{1}(iC1) & t_data.DrugName2=='-');

    for iC2 = 1:length(logC2)
        if iC1==1
            GrowthMx2(iC2) = t_data.RelGrowth(t_data.DrugName==Drugs(2) & ...
                t_data.Conc==Conc{2}(iC2) & t_data.DrugName2=='-');
        end

        val = t_data.RelGrowth(t_data.DrugName==Drugs(1) & ...
            t_data.Conc==Conc{1}(iC1) & t_data.DrugName2==Drugs(2) & ...
            t_data.Conc2==Conc{2}(iC2) );
        if ~isempty(val)
            GrowthMx(iC1,iC2) = val;
        end


    end
end
%%%%%%
% the combo should be an interpolation (like pcolor) in the log 10
% (for cases where the doses are not equidistant in the log domain)
%%%%%%


xlims = [1.5*logC2(1)-.6*logC2(2) 1.5*logC2(end)-.4*logC2(end-1)];
ylims = [1.5*logC1(1)-.6*logC1(2) 1.5*logC1(end)-.4*logC1(end-1)];


% small corner (control)
get_newaxes([pos(4,1) pos(4,2) xwidth ywidth])
imagesc(1,1,1,[-1 1])
set(gca,'xtick',1,'xticklabel',0, 'ytick',1,'yticklabel',0 , ...
    Plotting_parameters.axes{:},'ydir','normal')
xlim([.9 1.1])
ylim([.9 1.1])

% plot for Drug 1
get_newaxes(subpos(1,:),1)
imagesc(ones(1,sum((~isnan(GrowthMx1)))), logC1(~isnan(GrowthMx1)), ...
    GrowthMx1(~isnan(GrowthMx1)),[-1 1]);
ylabel([char(Drugs(1))  '[uM]'], Plotting_parameters.axislabel{:})
set(gca,'xtick',[], 'ytick', logC1(2:2:end), ...
    'yticklabel', num2cellstr(Conc{1}(2:2:end),'%.2g'), ...
    Plotting_parameters.axes{:},'ydir','normal')
% plot the GIxx for Drug 1
for i=1:length(GIvalues)
    plot([.5 1.5], log10(GIs(i,1))*[1 1],'-','color', GIcolors(i,:), ...
        'linewidth',2);
end
ylim(ylims)
xlim([.9 1.1])

% plot for Drug 2
get_newaxes(subpos(2,:),1)
imagesc(logC2(~isnan(GrowthMx2)), ones(1,sum(~isnan(GrowthMx2))), ...
    GrowthMx2(~isnan(GrowthMx2)),[-1 1]);
xlabel([char(Drugs(2)) '[uM]'], Plotting_parameters.axislabel{:})
set(gca,'ytick',[], 'xtick', logC2(2:2:end), ...
    'xticklabel', num2cellstr(Conc{2}(2:2:end),'%.2g'),...
    Plotting_parameters.axes{:},'ydir','normal')
% plot the GIxx for Drug 2
for i=1:length(GIvalues)
    plot(log10(GIs(i,2))*[1 1], [.5 1.5], '-','color', GIcolors(i,:), ...
        'linewidth',2);
end
xlim(xlims)
ylim([.9 1.1])

% plot for combination
get_newaxes(subpos(3,:),1)

% imagesc(logC2, logC1, GrowthMx,[-1 1]);
% interpolation
[idx1, idx2] = find(~isnan(GrowthMx));
intrplogC1 = (min(logC1(idx1))-.15):.15:(max(logC1(idx1))+.15);
intrplogC2 = (min(logC2(idx2))-.15):.15:(max(logC2(idx2))+.15);
% intrplogC1 = logC1(min(idx1):max(idx1));
% intrplogC2 = logC2(min(idx2):max(idx2));


intrpGrowth = scatteredInterpolant(logC1(idx1), logC2(idx2), ...
    GrowthMx(sub2ind(size(GrowthMx), idx1, idx2)));
[logC1pts, logC2pts] = meshgrid(intrplogC1, intrplogC2);

intrpGrowthMx = intrpGrowth(logC1pts(:), logC2pts(:));
intrpGrowthMx = reshape(intrpGrowthMx', length(intrplogC2), length(intrplogC1))';
imagesc(intrplogC2, intrplogC1, intrpGrowthMx, [-1 1]);

set(gca,'xtick',[], 'ytick',[],'ydir','normal')

% plot the ratio directions used for fit and the Aligned diretions for the fit
if ~isempty(t_rfits)
    plot([-5 5], ones(size(t_rfits.BAratio))*[-5 5]-t_rfits.BAratio*[1 1], ...
        '-', 'color', [.8 .8 .8]);
end
if ~isempty(AlignConc{1})
    plot([-5 5], log10(AlignConc{1})*[1 1], '--', 'color', [.8 .8 .8]);
end
if ~isempty(AlignConc{2})
    plot(log10(AlignConc{2})*[1 1], [-5 5], '--', 'color', [.8 .8 .8]);
end

% plot the isoGI cruves
h = zeros(length(GIvalues)+1,1);
for j=[ToRow(find(GIvalues~=50)) GI50idx]
    % theoric curve (additive)
    plot(log10(allGIs(j).NullConc(:,2)), log10(allGIs(j).NullConc(:,1)), ...
        '--','color', GIcolors(j,:), 'linewidth',1);
    % real curve
    plot(log10(allGIs(j).RealConc(:,2)), log10(allGIs(j).RealConc(:,1)), ...
        '.','color', GIcolors(j,:), 'markersize', 12);
    h(j) = plot(log10(allGIs(j).intrpConc(:,2)), log10(allGIs(j).intrpConc(:,1)), ...
        '-','color', GIcolors(j,:), 'linewidth',1);

end
h(end) = plot([NaN NaN], [NaN NaN], '--', 'color', [.3 .3 .3]);


% plot the ratio directions (1:1) and (GI50:GI50)
plot([min([logC1;logC2])-5 max([logC1;logC2])+5], ...
    [min([logC1;logC2])-5 max([logC1;logC2])+5], '-', 'color', [.5 .5 .5]);

plot([min([logC1;logC2])-5 max([logC1;logC2])+5]-log10(NormFactor(1)), ...
    [min([logC1;logC2])-5 max([logC1;logC2])+5]-log10(NormFactor(2)), '-k');


xlim(xlims)
ylim(ylims)
colormap(Plotting_parameters.cmapGrWBr)

hc = colorbar('location','east');
set(hc,'position',subpos(4,:))
ylabel(hc,'RelGrowth', Plotting_parameters.axislabel{:})
ylim(hc,[-.3 1])

hl = legend(h, [strcat('GI', num2cellstr(GIvalues)) 'CI=1'], ...
    Plotting_parameters.axislabel{:});
posl = get(hl,'position');
set(hl,'position', [subpos(5,1) subpos(5,2) posl(3:4)])
