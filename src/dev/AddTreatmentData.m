function t_annotated = AddTreatmentData(t_data, folder, fields)
%
%
%   annotate the data using the treatment files
%   adding the fields:
%       DrugName    (assuming only one drug per well)
%       Conc
%       any additional fields as given in input varaible 'fields'
%
%   variable folder is to specify where the Treatmentfiles are stored.
%



if ~exist('folder','var')
    folder = '';
end

%%
Trtfiles = setdiff(cellstr(unique(t_data.Treatmentfile)),'-');

% check for the treatment files
for iTf = 1:length(Trtfiles)
    if isempty(strfind(Trtfiles{iTf}, '.mat')), continue, end
    design_vars = load(fullfile(folder, Trtfiles{iTf}));
    designs = design_vars.design;
    for iD=1:size(designs)

    % avoid two drugs in the same well
    Drugs = reshape([designs{iD}.layout], [size(designs{iD}(1).layout) ...
        length(designs{iD})]);
    assert(all(all(sum(Drugs>0,3)<2)), 'some wells have two drugs')
    end
end

Conc = NaN(height(t_data),1);
DrugName = repmat({''}, height(t_data),1);

if exist('fields','var')
   datafields = cell(1, length(fields));
end


%%% this is wacky; need better
if exist('fields','var')
    warning('need better implementation of the design file')

    design_vars = load(fullfile(folder, Trtfiles{1}));
    designs = design_vars.design;
    idx = t_data.Treatmentfile=='-';

    LayoutIdx = sub2ind(size(designs{iD}(1).layout), t_data.Row(idx), t_data.Column(idx));

    for iF = 1:length(fields)
        temp = designs{1}(1).(fields{iF})(LayoutIdx);
        datafields{iF}(idx) = temp;
    end
end

for iTf = 1:length(Trtfiles)

    if ~isempty(strfind(Trtfiles{iTf}, '.tsv')) || ...
            ~isempty(strfind(Trtfiles{iTf}, '.txt'))
        % case of a .tsv file
        idx = t_data.Treatmentfile==Trtfiles{iTf};

        t_trt = TableToCategorical(tsv2table(fullfile(folder, Trtfiles{iTf})),0);

        [temp, idx2] = outerjoin(t_data(idx,:), t_trt, 'keys', {'Well'}, 'rightvariables', ...
            intersect(varnames(t_trt), {'HMSLid' 'DrugName' 'Conc'}), 'Type','left');
        temp = temp(sortidx(idx2),:);

        DrugName(idx) = cellstr(temp.DrugName);
        Conc(idx) = temp.Conc;

    elseif ~isempty(strfind(Trtfiles{iTf}, '.mat'))
        % case of a MATLAB file with a structure 'design'

        design_vars = load(fullfile(folder, Trtfiles{iTf}));
        designs = design_vars.design;

        DesignNumbers = setdiff(unique(t_data.DesignNumber),0);
        DesignNumbers = DesignNumbers(~isnan(DesignNumbers));
        for iD = DesignNumbers'
            idx = t_data.Treatmentfile==Trtfiles{iTf} & t_data.DesignNumber==iD;

            nDrugs = length(designs{iD});
            DrugNames = [{''}; ReducName({designs{iD}.name}, ',.')'];
            DrugConc = reshape([designs{iD}.layout], ...
                [size(designs{iD}(1).layout) nDrugs]);
            DrugIdx = sum((DrugConc>0).* ...
                repmat(reshape(1:nDrugs,1,1,[]),[size(designs{iD}(1).layout), 1]),3);
            DrugConc = sum(DrugConc,3);
            LayoutIdx = sub2ind(size(designs{iD}(1).layout), t_data.Row(idx), t_data.Column(idx));

            DrugName(idx) = DrugNames(DrugIdx(LayoutIdx)+1);
            Conc(idx) = DrugConc(LayoutIdx);

            if exist('fields','var')
                for iF = 1:length(fields)
                    temp = designs{iD}(1).(fields{iF})(LayoutIdx);
                    datafields{iF}(idx) = temp;
                end
            end

        end
    end
end
assert(all(strcmp(DrugName(Conc==0),'') | strcmp(DrugName(Conc==0),'DMSO')))

pert_type = repmat({''}, height(t_data),1);
pert_type(Conc==0 & ismember(t_data.Treatmentfile,Trtfiles)) = {'ctl_vehicle'};
pert_type(Conc>0 & ismember(t_data.Treatmentfile,Trtfiles)) = {'trt_cp'};


t_annotated = [t_data table(DrugName, Conc, pert_type) ];
if exist('fields','var')
    datafields = cellfun2(@ToColumn, datafields);
    t_annotated = [t_annotated table(datafields{:}, 'variablenames', fields)];
end

t_annotated = TableToCategorical(t_annotated,0);
