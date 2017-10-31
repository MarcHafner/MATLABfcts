function UpdatedClassifiedCells = CellCycleReClass(ClassifiedCells, SingleCelldata, ManualIdxs, SelectedWells, varargin)
% ClassifiedCells = CellCycleClassification(t_SingleCelldata, SingleCelldata, ...)
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
% output is a structure with:
%  t_results       -> table with all calculated data and metadata
%  allCellIdentity -> cell cycle phase assigned to each cell
%                       (dead are -1; unclassified 0; G1 1; S 2; G2 3; M 4)
%  t_qc            -> table with summary of qc results by well
%  t_summary       -> table with summary of qc results by condition
%  plotResults     -> structure with the values and gates for plotting EdU/DNA
%
%% assign inputs and prepare tests
assert(length(SingleCelldata)==height(ClassifiedCells.t_results))

if islogical(ManualIdxs)
    assert(length(SingleCelldata)==length(ManualIdxs))
    ManualIdxs = find(ManualIdxs);
    
elseif isvector(ManualIdxs)
    assert(all(ManualIdxs>0 & ManualIdxs<=length(SingleCelldata)))
    
else
    SelectedBarcodes = ManualIdxs;
    if ischar('SelectedBarcodes'), SelectedBarcodes = {SelectedBarcodes}; end
    if ischar('SelectedWells'), SelectedWells = {SelectedWells}; end
    assert(length(SelectedWells)==length(SelectedWells))
    
    ManualIdxs = false(length(SingleCelldata),1);
    for iB = 1:length(SelectedBarcodes)
        for i = 1:length(SelectedWells{iB})
            ManualIdxs(ClassifiedCells.t_results.Barcode==SelectedBarcodes{iB} & ...
                ClassifiedCells.t_results.Well==SelectedWells{i}) = true;
        end
    end
    ManualIdxs = find(ManualIdxs);
end

p = parseCellCycleInputs(varargin{:});

% flag for the modular approach
useEdU = ~isempty(p.Channelnames.EdU);
usepH3 = ~isempty(p.Channelnames.pH3);

% folder for storing data
if ~isempty(p.savefolder)
    grp_savefolder = [p.savefolder '/manual/'];
    mkdir(grp_savefolder)
else
    savefig = ''; % this will bypass saving
end

% unwrap the structure
t_results = ClassifiedCells.t_results;
t_qc = ClassifiedCells.t_qc;
allCellIdentity = ClassifiedCells.allCellIdentity;
plotResults = ClassifiedCells.plotResults;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop through the selected conditions

fprintf('\tTreatment wells (%i): ', length(ManualIdxs));
for iW = 1:length(ManualIdxs)
    
    if ismember(t_results.pert_type(ManualIdxs(iW)), p.NegCtrlLabel)
        welltag = 'Negctl_';
    elseif ismember(t_results.pert_type(ManualIdxs(iW)), p.PosCtrlLabel)
        welltag = 'Posctl_';
    else
        welltag = 'Trt_';
    end
    fprintf('% 2i/%i: %s - %s: %s, c=%.2g (%s)', iW, length(ManualIdxs), ...
        t_results.Barcode(ManualIdxs(iW)), t_results.Well(ManualIdxs(iW)), ...
        t_results.DrugName(ManualIdxs(iW)), t_results.Conc(ManualIdxs(iW)), welltag(1:end-1));
    logtxt = 'MANUAL:';
    
    
    LDRtxt = SingleCelldata(ManualIdxs(iW)).(p.Channelnames.LDR);
    DNA = SingleCelldata(ManualIdxs(iW)).(p.Channelnames.DNA);
    
    if isempty(LDRtxt), continue, end % display warning or print in log? <-----------------
    
    if ~isempty(p.savefolder), savefig = [grp_savefolder ...
            welltag char(t_results.Well(ManualIdxs(iW))) '_LDR_manual.jpg'];
    end
    [LiveCells, DeadCells, ~, CellOutcome, ~, ~, ltxt] = DeadCellFilter(LDRtxt, DNA, ...
        'savefig', savefig, 'interactive', true);
    % store the results
    logtxt = [logtxt ' ' ltxt];
    allCellIdentity{ManualIdxs(iW)} = CellOutcome;
    t_results(ManualIdxs(iW), {'LiveCells' 'DeadCells'}) = {LiveCells, DeadCells};
    % check and report outcome (special for toxic positive control)
    t_qc.PassFracDead(ManualIdxs(iW)) = true;
    
    if useEdU
        EdU = SingleCelldata(ManualIdxs(iW)).(p.Channelnames.EdU);
        if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                welltag char(t_results.Well(ManualIdxs(iW))) '_EdU_manual.jpg'];
        end
        [CCPks, CCfrac, ...
            plotResults(ManualIdxs(iW)).DNAGates, plotResults(ManualIdxs(iW)).EdUGates, ...
            CellIdentity, ...
            plotResults(ManualIdxs(iW)).logDNA, plotResults(ManualIdxs(iW)).logEdU,~,~,ltxt] = ...
            CCphases(DNA(CellOutcome==1), EdU(CellOutcome==1), ...
            'savefig', savefig, 'interactive', true);
        logtxt = [logtxt '; ' ltxt];
        
        if usepH3
            pH3 = SingleCelldata(ManualIdxs(iW)).(p.Channelnames.pH3);
            if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                    welltag char(t_results.Well(ManualIdxs(iW))) '_pH3_manual.jpg'];
            end
            [CCfrac, CellIdentity, ~, ~, ltxt] = pH3Filter(pH3(CellOutcome==1), ...
                CellIdentity, 'savefig', savefig, 'interactive', true);
            logtxt = [logtxt '; ' ltxt];
        end
        
        % store the results
        allCellIdentity{ManualIdxs(iW)}(CellOutcome==1) = CellIdentity;
        t_results(ManualIdxs(iW), {'CCfrac' 'CCPks'}) = {CCfrac {CCPks}};
        % check consistency and report outcome
        PassUnclass = mean(CellIdentity==0) < p.TestCutoffs.Unclass;
        PassCCphase = true;
        PassPksConsist = true;
        
        t_qc(ManualIdxs(iW), ...
            {'PassUnclass' 'PassCCphase' 'PassPksConsist'}) = ...
            {PassUnclass PassCCphase PassPksConsist};
        t_qc.log{ManualIdxs(iW)} = logtxt;
    end
    
end
fprintf(' done\n')



%%
t_summary = collapse(t_qc, @mean, ...
    'keyvars', unique([{'Barcode', 'CellLine'} p.ConditionKeys 'pert_type'],'stable'), ...
    'valvars', intersect(varnames(t_qc), {'PassFracDead' 'PassDeadConsist' 'PassUnclass' ...
    'PassCCphase' 'PassPksConsist' 'PassCCConsist'}, 'stable'));


UpdatedClassifiedCells = struct('t_results', {t_results}, 'allCellIdentity', {allCellIdentity}, ...
    't_qc', {t_qc}, 't_summary', {t_summary}, 'plotResults', {plotResults});