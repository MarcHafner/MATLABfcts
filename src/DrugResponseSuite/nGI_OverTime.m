function [t_nGITime, t_fitsTime] = nGI_OverTime(t_data, plate_keys, cond_keys, varargin)
% [t_nGITime, t_fitsTime] = nGI_OverTime(t_data, keys, varargin)
%   Normalized relative growth for different time intervals.
%   Sigmoidal fit on the drug response data (expect concentration in uM) to
%   extract the following parameters:
%       - GI50
%       - GIinf
%   All the outputs are saved in a table with annotations including drug
%   concentrations and initial values.
%
%   keys are used for aggregation of the data; default are : CellLine,
%   DrugName, Time, SeedingNumber, Date.
%
%   varargin:   - 'MinNdiv'     [0.1]
%               - 'MinDt'       [8 h]
%               - 'MaxDt'       
%               - 'minT0'       [0 h]
%               - 'T0date'      Input alternative to T0shift for timecourse:
%                                   date and time of the treatment
%               - 'pcutoff'     [0.1] cutoff for the p-value of a F-test against a flat line.
%               - 'forcefit'    [false], force a fit by adding value a low
%                                   concentrations,
%


p = inputParser;
addParameter(p, 'MinNDiv', 1/10, @isscalar);
addParameter(p, 'MinDt',   8, @isscalar);
addParameter(p, 'MaxDt',   96, @isscalar);
addParameter(p, 'minT0',    0, @isscalar);
addParameter(p, 'pcutoff', .1, @isscalar);
addParameter(p, 'forcefit', false, @isscalar);
parse(p,varargin{:})
p = p.Results;


if exist('plate_keys','var') && ~isempty(plate_keys)
    plate_keys = intersect(t_data.Properties.VariableNames, ...
        [{'CellLine' 'DrugName' 'Time' 'Date' 'Barcode' 'SeedingDensity'} plate_keys]);
else
    plate_keys = intersect(t_data.Properties.VariableNames, ...
        {'CellLine' 'DrugName' 'Time' 'Date' 'Barcode' 'SeedingDensity'});
end
if exist('cond_keys','var') && ~isempty(cond_keys)
    cond_keys = intersect(t_data.Properties.VariableNames, cond_keys);
else
    cond_keys = {};
end
    
t_keys = unique(t_data(t_data.DrugName~='-',setdiff([plate_keys cond_keys], {'Time' 'Date'})));


%%


t_fitsTime = table;
t_nGITime = table;
for ik = 1:height(t_keys)
    fprintf([strjoin(table2cellstr(t_keys(ik,:),0),'|') ' :']);
    %%
    subt = t_data(eqtable(t_keys(ik,:), t_data(:,[plate_keys cond_keys])),:);
    t_ctrl = sortrows(collapse(t_data(eqtable(t_keys(ik,:), t_data(:,setdiff(plate_keys,'DrugName'))) & ...
        t_data.pert_type=='ctl_vehicle' ,:), @mean, 'keyvars', plate_keys), 'Time');

    Times = t_ctrl.Time;
    assert(all(Times==unique(Times)));

    for iT = find(Times'>p.minT0)
        assert(t_ctrl.Time(iT)==Times(iT));
        % control
        Ctrl_DeltaT = (t_ctrl.Time((iT+1):end) - t_ctrl.Time(iT))/24;
        NDiv = log2(t_ctrl.Cellcount((iT+1):end)/t_ctrl.Cellcount(iT));
        Ctrl_AvDivRate = log2(t_ctrl.Cellcount((iT+1):end)/t_ctrl.Cellcount(iT))./Ctrl_DeltaT;

        idxEnd = find(Ctrl_DeltaT>=p.MinDt/24 & NDiv>=p.MinNDiv & Ctrl_DeltaT<=p.MaxDt/24);
        NDiv = NDiv(idxEnd);
        Ctrl_AvDivRate = Ctrl_AvDivRate(idxEnd);
        idxEnd = iT + idxEnd;
        fprintf(sprintf(' %.0f(%i);', Times(iT), length(idxEnd)));
        if mod(iT,10)==0, fprintf('\n'); end
        for iTE = 1:length(idxEnd)
            % treatment
            Conc = intersect(subt.Conc(subt.Time==Times(iT)), ...
                subt.Conc(subt.Time==Times(idxEnd(iTE))));
            Conc = setdiff(Conc, 0);
            t_trt = sortrows(collapse(subt(ismember(subt.Conc, Conc) & ...
                ismember(subt.Time, Times([iT idxEnd(iTE)])),:), ...
                @mean, 'keyvars', [plate_keys cond_keys 'Conc']), 'Time');

            nGI = NaN(length(Conc),1);
            parfor iC = 1:length(Conc)
                idx0 = t_trt.Time==Times(iT) & t_trt.Conc==Conc(iC);
                idxE = t_trt.Time==Times(idxEnd(iTE)) & t_trt.Conc==Conc(iC);
                trt_AvDivRate = log2(t_trt.Cellcount(idxE)/t_trt.Cellcount(idx0))/ ...
                    ((Times(idxEnd(iTE)) - Times(iT))/24);
                nGI(iC) = 2^(trt_AvDivRate/Ctrl_AvDivRate(iTE)) -1;
            end

            t_nGITime = [t_nGITime;
                [repmat([t_keys(ik,:), table(Times(iT), Times(idxEnd(iTE)),...
                diff(Times([iT idxEnd(iTE)])),  NDiv(iTE), ...
                'variablenames', {'T0' 'Tend' 'DeltaT' 'Ndiv'})], length(Conc),1) ...
                table(Conc, nGI)]];
            if nargout>1 && (length(Conc) > 4 || p.forcefit>0)
                if length(Conc) <= 4
                    n = 4-length(Conc)+ceil(p.forcefit);
                    Conc = [min(Conc).*(10.^(-n:-1)'); Conc];
                    nGI = [ones(n,1)*mean([1 max(nGI)]); nGI];
                end

                fitopt.pcutoff = p.pcutoff;
                [nGI50, ~, nGIinf, nGImax, nGIArea, nGI_r2, ~, nGI_fit] = ...
                    ICcurve_fit(Conc, nGI, 'nGI50', fitopt);
                t_fitsTime = [t_fitsTime;
                    [t_keys(ik,:) table(Times(iT), Times(idxEnd(iTE)), diff(Times([iT idxEnd(iTE)])), ...
                    NDiv(iTE), 'variablenames', {'T0' 'Tend' 'DeltaT' 'Ndiv'}), ...
                    table(nGI50, nGIinf, nGImax, nGIArea, nGI_r2) ...
                    table({nGI_fit}, {nGI'}, 'VariableNames', {'nGI_fit' 'nRelGrowth'})]];
            end
        end

    end
    fprintf('\n');
end

% matching the times to avoid rounding issues
uDt = unique(t_nGITime.DeltaT);
matchDt = [uDt cumsum([0;diff(uDt)>.5])];
for i=1:max(matchDt(:,2))
    t_nGITime.DeltaT(ismember(t_nGITime.DeltaT, matchDt(matchDt(:,2)==i,1))) = ...
        round(mean(matchDt(matchDt(:,2)==i,1)),2);
end
t_nGITime.Time = t_nGITime.T0 + t_nGITime.DeltaT/2;

if ~isempty(t_fitsTime)
    uDt = unique(t_fitsTime.DeltaT);
    matchDt = [uDt cumsum([0;diff(uDt)>.5])];
    for i=1:max(matchDt(:,2))
        t_fitsTime.DeltaT(ismember(t_fitsTime.DeltaT, matchDt(matchDt(:,2)==i,1))) = ...
            round(mean(matchDt(matchDt(:,2)==i,1)),2);
    end
    t_fitsTime.Time = t_fitsTime.T0 + t_fitsTime.DeltaT/2;
end
