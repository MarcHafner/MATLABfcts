function t_annotated = Annotate_CellCountData(t_data, folder, fields)
% t_annotated = Annotate_CellCountData(t_data, folder, fields)
%
%   annotate the data using the treatment files
%   adding the fields:
%       - DrugName    (assuming only one drug per well)
%       - Conc
%       - selected additional fields as given in input varaible 'fields' (by
%           default ALL 'Perturbations' in the array of structures Design)
%
%   variable folder is to specify where the TreatmentFiles are stored.
%


if ~exist('folder','var') || isempty(folder)
    folder = '';
end

%% import all treatment files and read the designs
Trtfiles = setdiff(cellstr(unique(t_data.TreatmentFile)),'-');

Designs = cell(length(Trtfiles),1);
correct_barcode = Designs;
for iTf = 1:length(Trtfiles)
    assert(exist(fullfile(folder, Trtfiles{iTf}), 'file')>0, ...
        'Treatment file %s missing in folder %s', Trtfiles{iTf}, folder)
    [~,~,ext] = fileparts(Trtfiles{iTf});
    fprintf('\tLoading %s\n', Trtfiles{iTf})
    
    if strcmp(ext,'.mat')
        temp = load(fullfile(folder, Trtfiles{iTf}));
        Designs{iTf} = temp.Design;
    elseif strcmp(ext,'.hpdd')
        [Designs{iTf}, correct_barcode{iTf}] = hpdd_importer(fullfile(folder, Trtfiles{iTf}));
        % because of redundant plates, barcodes have to be reassigned
        % barcode found in the hpdd files will overwrite other barcode
    else
        try % assume a tsv file
            Designs{iTf} = TextDesignFile_importer(fullfile(folder, Trtfiles{iTf}));
        catch err
            warnprintf(['expecting a .mat, .hpdd or a formated ' ...
                'tab-separated file as TreatmentFile\n'])
            rethrow(err)
        end
    end
    
end

%% look up all the drugs and perturbations in the designs

Ndrugs = 1;
Perturbations = {};
DrugNames = {};
t_HMSLids = table;
allDesigns = [Designs{:}];
for iD=1:size(allDesigns)
    % check for multiple drugs in the same well
    DrugConc = reshape([allDesigns(iD).Drugs.layout], [allDesigns(iD).plate_dims ...
        length(allDesigns(iD).Drugs)]);
    DrugNames = unique([DrugNames {allDesigns(iD).Drugs.DrugName}],'stable');
    if isfield(allDesigns(iD).Drugs,'HMSLid')
        t_HMSLids = [ t_HMSLids; table({allDesigns(iD).Drugs.DrugName}', ...
            {allDesigns(iD).Drugs.HMSLid}','VariableNames', {'DrugName' 'HMSLid'})];
    end
    if any(any(sum(DrugConc>0,3)>Ndrugs))
        Ndrugs = max(max(sum(DrugConc>0,3)));
        warnprintf('some wells have %i drugs, additional columns in output', ...
            Ndrugs)
    end
    
    % store all possible perturbations
    if isfield(allDesigns(iD), 'Perturbations')
        Perturbations = unique([Perturbations {allDesigns(iD).Perturbations.Name}], 'stable');
    end
end
t_HMSLids = unique(t_HMSLids);


%% declare the variables


% this is not the optimal way of storing multiple drugs because of the
% hierarcy between DrugName and Conc as well as the redudancy and possible
% swapping between Drug1 and Drug2 ; it makes matching between condition hard


if exist('fields','var')
    assert(all(ismember(fields, Perturbations)), ...
        'Not all ''fields'' found as perturbations in the design files')
else
    fields = Perturbations;
end
datafields = cell(1, length(fields));


%%


t_annotated = table;

for iTf = 1:length(Trtfiles)
    
    DesignNumbers = unique(t_data.DesignNumber(t_data.TreatmentFile==Trtfiles{iTf}));
    fprintf('Design %s:\n', Trtfiles{iTf} )
    for iDN = 1:length(DesignNumbers)
        if ~isempty(correct_barcode{iTf})
            DNidx = correct_barcode{iTf}.DesignNumber(DesignNumbers(iDN));
            assert(all(t_data.Barcode(t_data.DesignNumber==DesignNumbers(iDN))== ...
                correct_barcode{iTf}.Barcode(DesignNumbers(iDN))), ...
                'Mismatch between barcodes in the hpdd file and the DesignNumber')
            fprintf('\thpdd exp %i -> design %i\n', DesignNumbers(iDN), DNidx);
            
        else
            DNidx = iDN;
        end
        t_design = DrugDesignToTable(Designs{iTf}(DNidx), fields, DrugNames);
        idx = find(t_data.TreatmentFile==Trtfiles{iTf} & t_data.DesignNumber==DesignNumbers(iDN));
        [temp, ia] = innerjoin(t_data(idx,:), t_design, 'keys', 'Well');
        if height(temp)<sum(idx)
            warnprintf('Some wells (%s) have not annotations in file %s --> ignored', ...
                strjoin(cellstr(unique(t_data.Well(idx(setdiff(1:length(idx),ia))))'),', '), ...
                Trtfiles{iTf})
        end
        t_annotated = [t_annotated; temp];
    end
    
end

% fill up the columns for the untreated plates before merging
Untrtidx = t_data.TreatmentFile=='-';
if any(Untrtidx)
    NDrugs = sum(strfindcell(varnames(t_annotated),'DrugName')==1);
    pert_type = repmat({'Untrt'}, sum(Untrtidx),1);
    
    newvars = setdiff(varnames(t_annotated), [varnames(t_data) {'pert_type'}], 'stable');
    othervars = intersect(varnames(t_annotated), varnames(t_data), 'stable');
    
    temp = table;
    for iD=1:NDrugs
        % add the HMSLid by default; it will be filtered afterwards
        temp = [temp cell2table(repmat({'-'}, sum(Untrtidx),2), 'VariableName', ...
            {sprintf('DrugName%i', iD(iD>1)), sprintf('HMSLid%i', iD(iD>1))}) ...
            table(zeros(sum(Untrtidx),1), 'VariableName', {sprintf('Conc%i', iD(iD>1))})];
    end
    
    addvars = setdiff(newvars, varnames(temp));
    for i=1:length(addvars)
        %%% need to treat the case of SeedingNumber
        if strcmp(newvars{i}, 'SeedingNumber')
            error('Broadcasting of seeding number not implemented')
        end
        if isnumeric(t_annotated.(addvars{i}))
            temp = [temp table(zeros(sum(Untrtidx),1), 'VariableName', addvars(i))];
        else
            temp = [temp cell2table(repmat({'-'}, sum(Untrtidx),1), ...
                'VariableName', addvars(i))];
        end
    end
    
    
    t_annotated = [
        [t_data(Untrtidx,othervars) temp(:,newvars) table( pert_type)]
        t_annotated(:,othervars) t_annotated(:,newvars)  t_annotated(:,'pert_type')];
end

warnassert(height(t_annotated)==height(t_data), 'table went from %i to %i rows; check labels', ...
    height(t_data), height(t_annotated))
