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
%                       columns:  {'Barcode' 'CellLine' 'TrtFile' 'DesignNumber' 'Replicate' 'DrugSet'};
%
%                       column value examples:
%                           Barcode:    MN20140509001
%                           CellLine:   MDAMB231
%                           TrtFile: Design_20140307_Exp1_HME.mat
%                           DesignNumber: index of design in TrtFile
%                           Replicate:  replicate (treatment)
%                           DrugSet:  'A', 'B', ... or 'Day0' or 'Untrt'
%
%             FilesLocation.OutputFile = 'Results_20140310';
%             [needed in using 'savefile']
%                       File name for saving the results
%
%  savefile = save results in a file (need FilesLocation.OutputFile)
%


function t_final = import_ColumbusData_design(FilesLocation, extrafields)
if ~exist('extrafields','var')
    extrafields = {};
elseif ischar(extrafields)
    extrafields = {extrafields};
end

fields = {'folder' 'data_file' 'Barcode_table'};

for i=1:length(fields)
    if isfield(FilesLocation, fields{i})
        eval([fields{i} '=FilesLocation.' fields{i} ';'])
    else
        error('Field %s needed for the import', fields{i})
    end
end
Barcode_table = TableToString(Barcode_table, {'Barcode' 'CellLine' 'TrtFile' 'DrugSet'});
assert(all(ismember({'Barcode' 'CellLine' 'TrtFile' 'DesignNumber' 'Replicate' 'DrugSet'}, ...
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
TrtFile = cell(height(t_raw),1);
DrugSet = cell(height(t_raw),1);
Replicate = zeros(height(t_raw),1);
DesignNumber = zeros(height(t_raw),1);
Untrt = zeros(height(t_raw),1);
Day0 = zeros(height(t_raw),1);

cnt = 0;
for iBC = 1:height(Barcode_table)

    idx = find(strcmp(t_raw.Barcode, Barcode_table.Barcode{iBC}));
    assert(length(unique(PlateID(idx)))==1, '%i plate IDs found', ...
        length(unique(PlateID(idx))))
    assert(isempty(intersect(unique(PlateID(idx)), unique(PlateID(~idx)))))

    CellLine(idx) = Barcode_table.CellLine(iBC);
    Barcodes(idx) = Barcode_table.Barcode(iBC);
    TrtFile(idx) = Barcode_table.TrtFile(iBC);
    Replicate(idx) = Barcode_table.Replicate(iBC);
    DrugSet(idx) = Barcode_table.DrugSet(iBC);

    if strcmp('Untrt', Barcode_table.TrtFile(iBC))
        Untrt(idx) = 1;
    elseif strcmp('Day0', Barcode_table.TrtFile(iBC))
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

t_data = [table(CellLine,TrtFile,DesignNumber,DrugSet,Replicate,Untrt,Day0) t_raw(Usedidx,{'Well' 'Row' 'Column' 'Nuclei_NumberOfObjects'})];
t_data.Properties.VariableNames{'Nuclei_NumberOfObjects'} = 'Cellcount';
t_data = sortrows(TableToCategorical(t_data,{'CellLine' 'Well' 'DrugSet'}),...
    {'CellLine' 'TrtFile' 'DesignNumber' 'Column' 'Row'});
%%


welllist = cell(0,1);
colrowlist = NaN(0,2);
for i=1:24
    welllist(end+[1:16]) = strcat(num2cell('A':'P'), num2str(i,'%i'))';
    colrowlist(end+[1:16],:) = [ [1:16]' repmat(i,16,1)];
end

% get all drugs and treatments

t_Trt_Design = unique(t_data(t_data.DesignNumber>0,{'TrtFile' 'DesignNumber'}));
DrugNames = {};
for iTD=1:height(t_Trt_Design)

    load([folder t_Trt_Design.TrtFile{iTD}])

    temp = cellfun2(@(x) {x.name}, design(t_Trt_Design.DesignNumber(iTD)));
    temp = temp{:};
    DrugNames = unique([DrugNames MATLABsafename(ReducName(temp))]);
end
%%



t_CL_Trt_Design = sortrows(unique(t_data(t_data.DesignNumber>0,{'CellLine' 'TrtFile' 'DesignNumber'})));

DrugDoses = zeros(height(t_data), length(DrugNames));
Ctrl = false(height(t_data), 1);

for ix = extrafields
    eval([ix{:} '=zeros(height(t_data),1);']);
end

for iCTD = 1:height(t_CL_Trt_Design)

    data_idx = t_data.CellLine==t_CL_Trt_Design.CellLine(iCTD) & ...
        strcmp(t_data.TrtFile,t_CL_Trt_Design.TrtFile{iCTD}) & ...
        t_data.DesignNumber==t_CL_Trt_Design.DesignNumber(iCTD);

    Layout_idx = sub2ind([16 24], t_data.Row(data_idx), t_data.Column(data_idx));

    assert( length(unique(t_data.Replicate(data_idx)))==1 ,'Multiple replicates with same design')

    load([folder t_CL_Trt_Design.TrtFile{iCTD}])

    temp = cellfun2(@(x) {x.name}, design(t_CL_Trt_Design.DesignNumber(iCTD)));
    temp = MATLABsafename(ReducName(temp{:}));
    assert(all(ismember(temp, DrugNames)))

    for iD=1:length(temp)
        DrugDoses(data_idx, strcmp(DrugNames, temp{iD})) = ...
            design{t_CL_Trt_Design.DesignNumber(iCTD)}(iD).layout(Layout_idx);
    end
    for ix = extrafields
        eval([ix{:} '(data_idx)=design{t_CL_Trt_Design.DesignNumber(iCTD)}(iD).' ...
            ix{:} '(Layout_idx);']);
    end

    tempCtrl = false(16,24);
    tempCtrl(Layout_idx) = all(DrugDoses(data_idx,:)==0,2);
    if sum(tempCtrl)>20
        tempCtrl(sub2ind([16 24], [1 1 16 16], [1 24 1 24])) = false;
    end
    Ctrl(data_idx) = tempCtrl(Layout_idx);

end
%%
tempstr = 'Ctrl';
for ix = extrafields
    tempstr = [tempstr ',' ix{:}];
end

t_final = [t_data(:,{'CellLine' 'DrugSet' 'Replicate' 'Untrt' 'Day0' 'Well' 'Row' 'Column'}) ...
    eval(['table(' tempstr ')']) array2table(DrugDoses, 'variablenames', DrugNames) t_data(:,{'Cellcount'})];
