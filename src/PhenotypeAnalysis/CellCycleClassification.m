function [t_results, allLiveIdx, allCellIdentity, t_summary] = CellCycleClassification(t_SingleCelldata, SingleCelldata, varargin)
% [t_results, allLiveIdx, allCellIdentity, t_summary] = CellCycleClassification(t_SingleCelldata, SingleCelldata, ...)
%
% inputs are :
%  t_SingleCelldata  -> metadata for the plates/wells
%  SingleCelldata    -> object-level data (array of structures with one field per
%                               channel, same order as t_SingleCelldata)
%
% parameter inputs are:
%  Channelnames -> name of the variables in SingleCelldata (structure with
%                       field LDR, DNA, and optional: EdU, pH3
%  savefolder   -> name to save images of the results
        %%%   option to save only if a test fails ?? <<<<<<<<<<<<<<<----------------------
%  ConditionKeys-> name of keys to split conditions
%  PosCtrlLabel -> labels for the positive controls
%  NegCtrlLabel -> labels for the negative controls
%
% outputs are:
%  t_results       -> table with all calculated data and test results
%  allLiveIdx      -> boolean for live cells
%  allCellIdentity -> cell cycle phase assigned to each cell (dead are NaN)
%  t_summary       -> table with summary of test results by condition
%

%% assign inputs and prepare tests

assert(height(t_SingleCelldata)==length(SingleCelldata))

p = inputParser;

% names used for each channel
addParameter(p, 'Channelnames', ...
    struct('LDR', 'NucleiSelected_LDRTXTSERSpot8Px', ...
    'DNA', 'NucleiSelected_DNAcontent', ...
    'EdU', 'NucleiSelected_EdUINT', ...
    'pH3', 'NucleiSelected_PH3INT'), ...
    @(x) isstruct(x) && isfield(x,'LDR') && isfield(x,'DNA'))
% default to save the results as images; empty if not saving
addParameter(p, 'savefolder', 'temp/', @ischar)
% additional keys to consider for ananlyzing the results
addParameter(p, 'ConditionKeys', {}, @cellstr)
% flags for the type of treatments
addParameter(p, 'NegCtrlLabel', {'ctl_vehicle' 'untrt'}, @cellstr)
addParameter(p, 'PosCtrlLabel', {'ctl_toxic' 'ctl_G1' 'ctl_S' 'ctl_G2'}, @cellstr)
% define cutoffs for tests
addParameter(p, 'TestCutoffs', ...
     struct(...
    'FracDead', .1, ...
    'DeadConsist', .02, ...
    'Unclass', .05, ...
    'CCphase', .01, ...
    'CCphase_pos', @(x)1.3+x+.1, ...
    'CCConsist', .02, ...
    'PksConsist', log10(1.2)), ...
    @(x) isstruct(x) && all(isfield(x, {'FracDead', 'DeadConsist', 'Unclass', ...
    'CCphase', 'CCphase_pos', 'CCConsist', 'PksConsist'})))
addParameter(p, 'interactive', false, @islogical) %%%   option ?? <<<<<<<<<<<<<<<----------------------

parse(p,varargin{:});
p = p.Results;

% default to save the results as images; empty if not saving
p.savefolder = p.savefolder;
if ~isempty(p.savefolder), mkdir(p.savefolder); end
%%%   option to save only if a test fails ?? <<<<<<<<<<<<<<<----------------------


% flag for the modular approach
useEdU = ~isempty(p.Channelnames.EdU);
usepH3 = ~isempty(p.Channelnames.pH3);

%% %%%%%%%%%%%%%%%%%%
% start the analysis

% get all different conditions
t_groups = unique(t_SingleCelldata(:,unique([{'Barcode', 'CellLine'} p.ConditionKeys],'stable')));

% start the result report
t_results = [t_SingleCelldata(:,unique([{'Barcode', 'CellLine'} p.ConditionKeys ...
    {'pert_type' 'Well'}],'stable')) ...
    array2table(NaN(height(t_SingleCelldata),2), 'variablenames', ...
    {'LiveCells' 'DeadCells'}) ...
    array2table(false(height(t_SingleCelldata),2),'variablenames', ...
    {'PassFracDead' 'PassDeadConsist'})];

if useEdU
    t_results = [t_results ...
        table(NaN(height(t_SingleCelldata),4+usepH3), ...
        repmat({NaN(3,2)}, height(t_SingleCelldata),1), ...
        'variablenames', {'CCfrac' 'CCPks'}) ...
        array2table(false(height(t_SingleCelldata),4),'variablenames', ...
        {'PassUnclass' 'PassCCphase' 'PassCCConsist' 'PassPksConsist'})];
end

% structure to store all results per well
allLiveIdx = cell(height(t_results),1);
allCellIdentity = cell(height(t_results),1);
CellIdentity = [];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop through the different conditions

for iGr = 1:height(t_groups)
    fprintf('\n%i of %i: %s\n-----------------------------------\n', ...
        iGr, height(t_groups), strjoin(table2cellstr(t_groups(iGr,:),0),' '));
    
    % folder for storing data
    if ~isempty(p.savefolder)
        grp_savefolder = [p.savefolder ...
            strjoin(table2cellstr(t_groups(iGr,:),0),'_') '/'];
        mkdir(grp_savefolder)
    else
        savefig = ''; % this will bypass plotting and saving
    end
    
    %% %%%%%% all controls first %%%%%%%%%%%
    fprintf('\tAll controls\n');
    % first get all controls to define the dead cell gating
    NegCtrlidx = find(eqtable(t_groups(iGr,:), t_SingleCelldata) & ...
        ismember(t_SingleCelldata.pert_type, p.NegCtrlLabel));
    
    PosCtrlidx = find(eqtable(t_groups(iGr,:), t_SingleCelldata) & ...
        ismember(t_SingleCelldata.pert_type, p.PosCtrlLabel));
    
    % if no positive control, use all wells            <<<<<<<<< -- not sure if right --------
    if any(PosCtrlidx), AllCtrlidx = [NegCtrlidx; PosCtrlidx]; else
        AllCtrlidx = find(eqtable(t_groups(iGr,:), t_SingleCelldata)); end
    %  <<<<<<<<< -- maybe safer to only use the negative controls --------
    
    LDRtxt = vertcat(SingleCelldata(AllCtrlidx).(p.Channelnames.LDR));
    DNA = vertcat(SingleCelldata(AllCtrlidx).(p.Channelnames.DNA));
    
    if ~isempty(p.savefolder), savefig = [grp_savefolder 'All_ctl_LDR.jpg']; end
    [~, ~, RefDeadGates, ~, LDRlims, DNAlims] = DeadCellFilter(LDRtxt, DNA, ...
        'savefig', savefig);
    
    
    % run on all negative controls
    LDRtxt = vertcat(SingleCelldata(NegCtrlidx).(p.Channelnames.LDR));
    DNA = vertcat(SingleCelldata(NegCtrlidx).(p.Channelnames.DNA));
    if ~isempty(p.savefolder), savefig = [grp_savefolder 'All_Negctl_LDR.jpg']; end
    [~, RefFracDead, ~, LiveIdx] = DeadCellFilter(LDRtxt, DNA, ...
        'savefig', savefig, 'LDRlims', LDRlims, 'DNAlims', DNAlims);
    RefFracDead = RefFracDead/length(LDRtxt);
    if useEdU
        EdU = vertcat(SingleCelldata(NegCtrlidx).(p.Channelnames.EdU));
        if ~isempty(p.savefolder), savefig = [grp_savefolder 'All_Negctl_EdU.jpg']; end
        [RefCCPks, RefCCfrac,~,~,CellIdentity,~,~,~, EdUlims] = CCphases(DNA(LiveIdx), ...
            EdU(LiveIdx), 'savefig', savefig, 'DNAlims', DNAlims);
    end
    if usepH3
        pH3 = vertcat(SingleCelldata(NegCtrlidx).(p.Channelnames.pH3));
        if ~isempty(p.savefolder), savefig = [grp_savefolder 'All_Negctl_pH3.jpg']; end
        [RefCCfrac, ~, RefpH3cutoff, pH3lims] = pH3Filter(pH3(LiveIdx), CellIdentity, ...
            'savefig', savefig);
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check how it applies for each negative control
    fprintf('\tNeg controls (%i): ', length(NegCtrlidx));
    for iW = 1:length(NegCtrlidx)
        fprintf(' %s', char(t_SingleCelldata.Well(NegCtrlidx(iW))));
        
        LDRtxt = SingleCelldata(NegCtrlidx(iW)).(p.Channelnames.LDR);
        DNA = SingleCelldata(NegCtrlidx(iW)).(p.Channelnames.DNA);
        
        if isempty(LDRtxt), continue, end % display warning or print in log? <-----------------
        
        if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                'Negctl_' char(t_SingleCelldata.Well(NegCtrlidx(iW))) '_LDR.jpg'];
        end
        [LiveCells, DeadCells, ~, LiveIdx] = DeadCellFilter(LDRtxt, DNA, ...
            'Gates', RefDeadGates, 'LDRlims', LDRlims, 'DNAlims', DNAlims, ...
            'savefig', savefig);
        % store the results
        allLiveIdx{NegCtrlidx(iW)} = find(LiveIdx);
        t_results(NegCtrlidx(iW), {'LiveCells' 'DeadCells'}) = {LiveCells, DeadCells};
        % check consistency and report outcome
        PassFracDead = DeadCells/length(LDRtxt) < p.TestCutoffs.FracDead;
        PassDeadConsist = abs(DeadCells/length(LDRtxt) - RefFracDead) < p.TestCutoffs.DeadConsist;
        
        t_results(NegCtrlidx(iW), {'PassFracDead' 'PassDeadConsist'}) = ...
            {PassFracDead PassDeadConsist};
        
        if useEdU
            EdU = SingleCelldata(NegCtrlidx(iW)).(p.Channelnames.EdU);
            if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                    'Negctl_' char(t_SingleCelldata.Well(NegCtrlidx(iW))) '_EdU.jpg'];
            end
            [CCPks, CCfrac,~,~,CellIdentity] = CCphases(DNA(LiveIdx), EdU(LiveIdx), ...
                'savefig', savefig, 'DNAlims', DNAlims, 'EdUlims', EdUlims);
            
            if usepH3
                pH3 = SingleCelldata(NegCtrlidx(iW)).(p.Channelnames.pH3);
                if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                        'Negctl_' char(t_SingleCelldata.Well(NegCtrlidx(iW))) '_pH3.jpg'];
                end
                CCfrac = pH3Filter(pH3(LiveIdx), ...
                    CellIdentity, 'savefig', savefig, 'pH3lims', pH3lims, ...
                    'pH3cutoff', RefpH3cutoff);
            end
            
            % store the results
            allCellIdentity{NegCtrlidx(iW)} = NaN(length(DNA),1);
            allCellIdentity{NegCtrlidx(iW)}(allLiveIdx{NegCtrlidx(iW)}) = CellIdentity;
            t_results(NegCtrlidx(iW), {'CCfrac' 'CCPks'}) = {CCfrac {CCPks}};
            % check consistency and report outcome
            PassUnclass = mean(CellIdentity==0) < p.TestCutoffs.Unclass;
            PassCCphase = all(CCfrac(1:3)>.05) & CCfrac(end-1) > p.TestCutoffs.CCphase;
            PassCCConsist = max(abs(CCfrac - RefCCfrac)) < p.TestCutoffs.CCConsist;
            PassPksConsist = max(abs(CCPks(:,1)-RefCCPks(:,1))) < p.TestCutoffs.PksConsist;
            
            t_results(NegCtrlidx(iW), ...
                {'PassUnclass' 'PassCCphase' 'PassCCConsist' 'PassPksConsist'}) = ...
                {PassUnclass PassCCphase PassCCConsist PassPksConsist};
        end
    end
    fprintf(' done\n')
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check how it applies for each positive control
    
    fprintf('\tPos controls (%i): ', length(PosCtrlidx));
    for iW = 1:length(PosCtrlidx)
        fprintf(' %s', char(t_SingleCelldata.Well(PosCtrlidx(iW))));
        
        LDRtxt = SingleCelldata(PosCtrlidx(iW)).(p.Channelnames.LDR);
        DNA = SingleCelldata(PosCtrlidx(iW)).(p.Channelnames.DNA);
        
        if isempty(LDRtxt), continue, end % display warning or print in log? <-----------------
        
        if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                'Posctl_' char(t_SingleCelldata.Well(PosCtrlidx(iW))) '_LDR.jpg'];
        end
        [LiveCells, DeadCells, ~, LiveIdx] = DeadCellFilter(LDRtxt, DNA, ...
            'Gates', RefDeadGates, 'LDRlims', LDRlims, 'DNAlims', DNAlims, ...
            'savefig', savefig);
        % store the results
        allLiveIdx{PosCtrlidx(iW)} = find(LiveIdx);
        t_results(PosCtrlidx(iW), {'LiveCells' 'DeadCells'}) = {LiveCells, DeadCells};
        % check and report outcome (special for toxic positive control)
        t_results.PassFracDead(PosCtrlidx(iW)) = (DeadCells/length(LDRtxt)) > ...
            (RefFracDead + .5*(t_SingleCelldata.pert_type(PosCtrlidx(iW))=='ctl_toxic'));
        
        if useEdU
            EdU = SingleCelldata(PosCtrlidx(iW)).(p.Channelnames.EdU);
            if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                    'Posctl_' char(t_SingleCelldata.Well(PosCtrlidx(iW))) '_EdU.jpg'];
            end
            [CCPks, CCfrac,~,~,CellIdentity] = CCphases(DNA(LiveIdx), EdU(LiveIdx), ...
                'savefig', savefig, 'DNAlims', DNAlims, 'EdUlims', EdUlims);
            
            if usepH3
                pH3 = SingleCelldata(PosCtrlidx(iW)).(p.Channelnames.pH3);
                if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                        'Posctl_' char(t_SingleCelldata.Well(PosCtrlidx(iW))) '_pH3.jpg'];
                end
                CCfrac = pH3Filter(pH3(LiveIdx), ...
                    CellIdentity, 'savefig', savefig, 'pH3lims', pH3lims, ...
                    'pH3cutoff', RefpH3cutoff);
            end
            
            % store the results
            allCellIdentity{PosCtrlidx(iW)} = NaN(length(DNA),1);
            allCellIdentity{PosCtrlidx(iW)}(allLiveIdx{PosCtrlidx(iW)}) = CellIdentity;
            t_results(PosCtrlidx(iW), {'CCfrac' 'CCPks'}) = {CCfrac CCPks};
            % check consistency and report outcome
            PassUnclass = mean(CellIdentity==0) < .05;
            switch t_SingleCelldata.pert_type(PosCtrlidx(iW))
                case ctl_G1
                    PassCCphase = CCfrac(1) > p.TestCutoffs.CCphase_pos(RefCCfrac(1));
                case ctl_S
                    PassCCphase = CCfrac(2) > p.TestCutoffs.CCphase_pos(RefCCfrac(2));
                case ctl_G2
                    PassCCphase = CCfrac(3) > p.TestCutoffs.CCphase_pos(RefCCfrac(3));
                otherwise
                    PassCCphase = true;  % <<<<<<<<<<<<------------ need to check
            end
            PassPksConsist = max(abs(CCPks(:,1)-RefCCPks(:,1))) < p.TestCutoffs.PksConsist;
            
            t_results(PosCtrlidx(iW), ...
                {'PassUnclass' 'PassCCphase' 'PassPksConsist'}) = ...
                {PassUnclass PassCCphase PassPksConsist};
        end
    end
    fprintf(' done\n')
    
    %% check consistency among the positive controls with the same treatments
    t_Posctrl = unique(t_SingleCelldata(PosCtrlidx, ...
        unique([{'Barcode', 'CellLine'} p.ConditionKeys {'DrugName' 'Conc'}])));
    for iP=1:height(t_Posctrl)
        idx = find(eqtable(t_Posctrl(iP,:), t_SingleCelldata));
        FracDead = t_results.DeadCells(idx)./ ...
            (t_results.DeadCells(idx) + t_results.LiveCells(idx));
        t_results.PassDeadConsist(idx) = abs(FracDead - mean(FracDead)) < p.TestCutoffs.DeadConsist;
        
        CCfrac = t_results.CCfrac(idx,:);
        t_results.PassCCConsist(idx) = max(abs(CCfrac - ...
            repmat(mean(CCfrac,1),length(idx),1)),[],2) < 2*p.TestCutoffs.CCConsist;
    end
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % work out all the other wells in the same condition
    Trtidx = find(eqtable(t_groups(iGr,:), t_SingleCelldata) & ...
        ~ismember(t_SingleCelldata.pert_type, [p.NegCtrlLabel p.PosCtrlLabel]));
    fprintf('\tTreatment wells (%i): ', length(Trtidx));
    for iW = 1:length(Trtidx)
        if mod(iW, ceil(length(Trtidx)/8))==0, fprintf(' %i', iW); end
        
        LDRtxt = SingleCelldata(Trtidx(iW)).(p.Channelnames.LDR);
        DNA = SingleCelldata(Trtidx(iW)).(p.Channelnames.DNA);
        
        if isempty(LDRtxt), continue, end % display warning or print in log? <-----------------
        
        if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                'Trt_' char(t_SingleCelldata.Well(Trtidx(iW))) '_LDR.jpg'];
        end
        [LiveCells, DeadCells, ~, LiveIdx] = DeadCellFilter(LDRtxt, DNA, ...
            'Gates', RefDeadGates, 'LDRlims', LDRlims, 'DNAlims', DNAlims, ...
            'savefig', savefig);
        % store the results
        allLiveIdx{Trtidx(iW)} = find(LiveIdx);
        t_results(Trtidx(iW), {'LiveCells' 'DeadCells'}) = {LiveCells, DeadCells};
        % check and report outcome (special for toxic positive control)
        t_results.PassFracDead(Trtidx(iW)) = true;
        
        if useEdU
            EdU = SingleCelldata(Trtidx(iW)).(p.Channelnames.EdU);
            if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                    'Trt_' char(t_SingleCelldata.Well(Trtidx(iW))) '_EdU.jpg'];
            end
            [CCPks, CCfrac,~,~,CellIdentity] = CCphases(DNA(LiveIdx), EdU(LiveIdx), ...
                'savefig', savefig, 'DNAlims', DNAlims, 'EdUlims', EdUlims);
            
            if usepH3
                pH3 = SingleCelldata(Trtidx(iW)).(p.Channelnames.pH3);
                if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                        'Trt_' char(t_SingleCelldata.Well(Trtidx(iW))) '_pH3.jpg'];
                end
                CCfrac = pH3Filter(pH3(LiveIdx), ...
                    CellIdentity, 'savefig', savefig, 'pH3lims', pH3lims, ...
                    'pH3cutoff', RefpH3cutoff);
            end
            
            % store the results
            allCellIdentity{Trtidx(iW)} = NaN(length(DNA),1);
            allCellIdentity{Trtidx(iW)}(allLiveIdx{Trtidx(iW)}) = CellIdentity;
            t_results(Trtidx(iW), {'CCfrac' 'CCPks'}) = {CCfrac {CCPks}};
            % check consistency and report outcome
            PassUnclass = mean(CellIdentity==0) < p.TestCutoffs.Unclass;
            PassCCphase = true;  
            PassPksConsist = max(abs(CCPks(:,1)-RefCCPks(:,1))) < p.TestCutoffs.PksConsist;
            
            t_results(Trtidx(iW), ...
                {'PassUnclass' 'PassCCphase' 'PassPksConsist'}) = ...
                {PassUnclass PassCCphase PassPksConsist};
        end
    end
    fprintf(' done\n')
    
    %% check consistency among the wells with same treatments
    t_trt = unique(t_SingleCelldata(Trtidx,  ...
        unique([{'Barcode', 'CellLine'} p.ConditionKeys {'DrugName' 'Conc'}])));
    for iP = 1:height(t_trt)
        idx = find(eqtable(t_trt(iP,:), t_SingleCelldata));
        FracDead = t_results.DeadCells(idx)./ ...
            (t_results.DeadCells(idx) + t_results.LiveCells(idx));
        t_results.PassDeadConsist(idx) = abs(FracDead - nanmean(FracDead)) < p.TestCutoffs.DeadConsist;
        
        CCfrac = t_results.CCfrac(idx,:);
        t_results.PassCCConsist(idx) = max(abs(CCfrac - ...
            repmat(nanmean(CCfrac,1),length(idx),1)),[],2) < 2*p.TestCutoffs.CCConsist;
    end
    
    fprintf('\tdone\n')
end


%%
t_summary = collapse(t_results, @mean, ...
    'keyvars', unique([{'Barcode', 'CellLine'} p.ConditionKeys 'pert_type'],'stable'), ...
    'valvars', intersect(varnames(t_results), {'PassFracDead' 'PassDeadConsist' 'PassUnclass' ...
    'PassCCphase' 'PassPksConsist' 'PassCCConsist'}, 'stable'));

