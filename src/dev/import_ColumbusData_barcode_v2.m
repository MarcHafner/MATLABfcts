% [t_CL, Drugs, CellLines] = import_ColumbusData_barcode(FilesLocation, savefile)
%
%   import the results of the D300 treated plates scanned with Columbus
%
%   all path to files and parameters have to be passed in FilesLocation as:
%             FilesLocation.folder = './results_20140307/';
%                       folder with the designs and the results (in onefile)
%             FilesLocation.data_file = '20140307_Exp1_Columbusoutput.txt';
%                       result file (output from Columbus)
%             FilesLocation.Barcode_table = table with the following
%                       columns:  Barcode   CellLine    Designfile  Replicate  Untrt   Day0
%
%                       column value examples:
%                           Barcode:    MN20140509001
%                           CellLine:   MDAMB231
%                           Designfile: Design_20140307_Exp1_HME.mat
%                                           (or '' if untreated)
%                           Replicate:  replicate (treatment)
%                           Untrt:  false / true
%                           Day0:  false / true
%
%             FilesLocation.OutputFile = 'Results_20140310';
%             [needed in using 'savefile']
%                       File name for saving the results
%
%  savefile = save results in a file (need FilesLocation.OutputFile)
%


function [t_CL, Drugs, CellLines] = import_ColumbusData_barcode(FilesLocation, savefile)

if ~exist('savefile','var'),
    savefile = 0;
end

fields = {'folder' 'data_file' 'Barcode_table' 'OutputFile'};
Mandatory = [1   1   1  savefile];

for i=1:length(fields)
    if isfield(FilesLocation, fields{i})
        eval([fields{i} '=FilesLocation.' fields{i} ';'])
    elseif Mandatory(i)
        error('Field %s needed for the import', fields{i})
    else
        eval([fields{i} '={};'])
    end
end
assert(all(ismember({'Barcode' 'CellLine' 'TrtFile' 'DesignNumber'}, ...
    Barcode_table.Properties.VariableNames)))

%%

t_raw = readtable([folder data_file],'delimiter','\t');

if length(unique(t_raw.NumberOfAnalyzedFields))>1
    Nref = median(t_raw.NumberOfAnalyzedFields);
    warning('wells with missing fields: %i', ...
        unique(t_raw.NumberOfAnalyzedFields))
    for i = find(t_raw.NumberOfAnalyzedFields~=Nref)'
        t_raw.Nuclei_NumberOfObjects(i) = t_raw.Nuclei_NumberOfObjects(i)*...
            (Nref/t_raw.NumberOfAnalyzedFields(i));
    end
end

Code_date = regexp(t_raw.Result,' > ','split');
Code_date = vertcat(Code_date{:});
t_raw = [cell2table(Code_date,'VariableName',{'Barcode' 'date'}) t_raw(:,2:end)];

assert(length(unique(t_raw.Barcode))>=height(Barcode_table), ...
    'Found %i plates, but expected %i plates', ...
    length(unique(t_raw.Barcode)), height(Barcode_table));

url = regexp(t_raw.URL,'/','split');
PlateID = vertcat(url{:});
PlateID = cellfun(@str2num,PlateID(:,end-1));

CellLine = cell(height(t_raw),1);
Barcodes = cell(height(t_raw),1);
Trt = cell(height(t_raw),1);
Repeat = cell(height(t_raw),1);
DesignNumber = zeros(height(t_raw),1);
Untrt = zeros(height(t_raw),1);
Day0 = zeros(height(t_raw),1);

cnt = 0;
for iBC = 1:height(Barcode_table)

    idx = find(strcmp(t_raw.Barcode, Barcode_table.Barcode{iBC}));
    assert(length(unique(PlateID(idx)))==1, '%i plate IDs found', ...
        length(unique(PlateID(idx))))
    assert(isempty(intersect(unique(PlateID(idx)), unique(PlateID(~idx)))))

    CellLine(idx) = {char(Barcode_table.CellLine(iBC))};
    Barcodes(idx) = Barcode_table.Barcode(iBC);
    Trt(idx) = {char(Barcode_table.TrtFile(iBC))};
    Repeat(idx) = {char(Barcode_table.Repeat(iBC))};

    if strcmp('Untrt', char(Barcode_table.TrtFile(iBC)))
        Untrt(idx) = 1;
    elseif strcmp('Day0', char(Barcode_table.TrtFile(iBC)))
        Day0(idx) = 1;
    else
        DesignNumber(idx) = Barcode_table.DesignNumber(iBC);
    end

    assert(all(DesignNumber(idx)>0 | Untrt(idx)>0 | Day0(idx)>0))
    cnt = cnt+1;

end
assert(cnt<=length(unique(PlateID)))
if cnt<length(unique(PlateID))
    warning('some entries in the result file are unused!')
    Usedidx = ~cell2mat(cellfun2(@isempty,CellLine));
    CellLine = CellLine(Usedidx);
    DesignNumber = DesignNumber(Usedidx);
    Untrt = Untrt(Usedidx);
    Day0 = Day0(Usedidx);
else
    Usedidx = 1:height(t_raw);
end

t_data = [table(CellLine,DesignNumber,Untrt,Day0) t_raw(Usedidx,{'Well' 'Row' 'Column' 'Nuclei_NumberOfObjects'})];
t_data.Properties.VariableNames{'Nuclei_NumberOfObjects'} = 'Cellcount';
t_data = TableToCategorical(t_data,[1 5]);
%%
CellLines = categories(t_data.CellLine);
t_CL = cell(1,length(CellLines));
Drugs = t_CL;

welllist = cell(0,1);
colrowlist = NaN(0,2);
for i=1:24
    welllist(end+[1:16]) = strcat(num2cell('A':'P'), num2str(i,'%i'))';
    colrowlist(end+[1:16],:) = [ [1:16]' repmat(i,16,1)];
end

for iCL=1:length(CellLines)

    t_CL{iCL} = t_data(t_data.CellLine==CellLines{iCL},:);

    for iDe = unique(t_CL{iCL}.DesignNumber)'
        not_measured = ~ismember(welllist, t_CL{iCL}.Well(t_CL{iCL}.DesignNumber==iDe));
        temp = cell2table([repmat(CellLines(iCL), sum(not_measured),1), ...
            repmat({iDe NaN NaN}, sum(not_measured),1) ...
            welllist(not_measured) num2cell(colrowlist(not_measured,:)) ...
            repmat({NaN}, sum(not_measured),1)], ...
            'variablenames', t_CL{iCL}.Properties.VariableNames);
        temp = TableToCategorical(temp,[1 5]);
        t_CL{iCL} = [t_CL{iCL};temp];
    end

    tempDrugs = {};
    IsPrimary = {};
    SingleDoses = {};
    ComboDoses = {};
    for iDe = find(Barcode_table.CellLine==CellLines(iCL) &...
            Barcode_table.DesignNumber>0)'

        load([folder char(Barcode_table.TrtFile(iDe))])

        temp = cellfun2(@(x) {x.name}, design(Barcode_table.DesignNumber(iDe)));
        temp = temp{:};
        for i=1:length(temp)
            if ~ismember(temp{i}, tempDrugs)
                tempDrugs =  [tempDrugs temp{i}];
                if isfield(design{Barcode_table.DesignNumber(iDe)},'Primary')
                    IsPrimary{end+1} = design{Barcode_table.DesignNumber(iDe)}.Primary;
                end
                SingleDoses{end+1} = design{Barcode_table.DesignNumber(iDe)}.SingleDoses;

                if isfield(design{Barcode_table.DesignNumber(iDe)},'Doses')
                    ComboDoses{end+1} = design{Barcode_table.DesignNumber(iDe)}.Doses;
                end
                if isfield(design{Barcode_table.DesignNumber(iDe)},'ComboDoses')
                    ComboDoses{end+1} = design{Barcode_table.DesignNumber(iDe)}.ComboDoses;
                end
            end
        end
    end
    if isempty(IsPrimary) && isempty(ComboDoses)
        Drugs{iCL} = struct('DrugName',MATLABsafename(ReducName(tempDrugs)), ...
            'SingleDoses',SingleDoses);
    elseif isempty(IsPrimary)
        Drugs{iCL} = struct('DrugName',MATLABsafename(ReducName(tempDrugs)), ...
            'SingleDoses',SingleDoses, ...
            'ComboDoses', ComboDoses);
    else
        Drugs{iCL} = struct('DrugName',MATLABsafename(ReducName(tempDrugs)), ...
            'Primary',IsPrimary, ...
            'SingleDoses',SingleDoses, ...
            'ComboDoses', ComboDoses);
    end


    DrugDoses = zeros(height(t_CL{iCL}),length(tempDrugs));
    t_CL{iCL} = sortrows(t_CL{iCL},{'DesignNumber' 'Column' 'Row'});

    for iDe = find(Barcode_table.CellLine==CellLines(iCL) &...
            Barcode_table.DesignNumber>0)'


        load([folder char(Barcode_table.TrtFile(iDe))])

        temp = cellfun2(@(x) {x.name}, design(Barcode_table.DesignNumber(iDe)));
        temp = temp{:};
        assert(all(ismember(temp, tempDrugs)))
        for iD=1:length(temp)
            idx = sub2ind([16 24], t_CL{iCL}.Row(t_CL{iCL}.DesignNumber==iDe & ~isnan(t_CL{iCL}.Cellcount)), ...
                t_CL{iCL}.Column(t_CL{iCL}.DesignNumber==iDe & ~isnan(t_CL{iCL}.Cellcount)) );
            DrugDoses(t_CL{iCL}.DesignNumber==iDe & ~isnan(t_CL{iCL}.Cellcount), ...
                strcmp(tempDrugs, design{Barcode_table.DesignNumber(iDe)}(iD).name)) = ...
                design{Barcode_table.DesignNumber(iDe)}(iD).layout(idx);
        end


    end


    CtrlIdx = all(DrugDoses==0,2);
    for iDe = setdiff(unique(t_CL{iCL}.DesignNumber), 0)'
        % remove corner as control to avoid bias (only the case if more
        % than 20 controls overall (enough control remain)
        CornerIdx = find(t_CL{iCL}.DesignNumber==iDe & ...
            ismember(t_CL{iCL}.Well,{'A1' 'A24' 'P1' 'P24'}));
        if all(CtrlIdx(CornerIdx)==1)
            CtrlIdx(CornerIdx)=0;
        end
    end

    t_CL{iCL} = [t_CL{iCL} array2table([DrugDoses CtrlIdx], ...
        'VariableNames',{Drugs{iCL}.DrugName 'Ctrl'})];

    t_CL{iCL} = sortrows(t_CL{iCL},{'DesignNumber' 'Column' 'Untrt' ...
        'Day0' 'Ctrl'});
end

if savefile
    save([OutputFile '.mat'], 't_CL', 'Drugs', 'CellLines')
end
