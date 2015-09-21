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
assert(all(ismember({'Barcode' 'CellLine' 'DesignFile'}, ...
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
DesignSet = cell(height(t_raw),1);
Replicate = zeros(height(t_raw),1);
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
    DesignSet(idx) = Barcode_table.DesignSet(iBC);
    Replicate(idx) = Barcode_table.Replicate(iBC);

    if ismember('Untrt', Barcode_table.Properties.VariableNames)
        Untrt(idx) = Barcode_table.Untrt(iBC);
    end

    if ismember('Day0',Barcode_table.Properties.VariableNames)
        Day0(idx) = Barcode_table.Day0(iBC);
    end

    assert(all(Replicate(idx)>0 | Untrt(idx)>0 | Day0(idx)>0))
    cnt = cnt+1;

end
assert(cnt<=length(unique(PlateID)))
if cnt<length(unique(PlateID))
    warning('some entries in the result file are unused!')
    Usedidx = ~cell2mat(cellfun2(@isempty,CellLine));
    CellLine = CellLine(Usedidx);
    Barcodes = Barcodes(Usedidx);
    DesignSet = DesignSet(Usedidx);
    Replicate = Replicate(Usedidx);
    Untrt = Untrt(Usedidx);
    Day0 = Day0(Usedidx);
else
    Usedidx = 1:height(t_raw);
end

CellLine = categorical(CellLine);
t_data = [table(CellLine,Barcodes,DesignSet,Replicate,Untrt,Day0) t_raw(Usedidx,{'Well' 'Row' 'Column' 'Nuclei_NumberOfObjects'})];
t_data.Properties.VariableNames{'Nuclei_NumberOfObjects'} = 'Cellcount';

%%
CellLines = categories(CellLine);
t_CL = cell(1,length(CellLines));
Drugs = t_CL;
DrugsNames = t_CL;

for iCL=1:length(CellLines)

    t_CL{iCL} = t_data(t_data.CellLine==CellLines{iCL},:);

    Designs = setdiff(unique(Barcode_table.DesignSet(...
        strcmp(Barcode_table.CellLine,CellLines{iCL}))),'');

    DrugsNames{iCL} = {};
    Drugs{iCL} = [];
    for iDe=1:length(Designs)
        design_file = Barcode_table.DesignFile(strcmp(Barcode_table.CellLine,...
            CellLines(iCL)) & Barcode_table.Replicate>0 & strcmp(Barcode_table.DesignSet, Designs{iDe}));
        load([folder design_file{1}])

        DrugsNames{iCL} = [ DrugsNames{iCL} ; {drugs_struct.name}'];

        if isfield(drugs_struct,'Primary')
            temp = struct('DrugName',MATLABsafename(ReplaceName({drugs_struct.name})), ...
                'Primary',{drugs_struct.Primary}, ...
                'SingleDoses',{drugs_struct.SingleDoses});
        else
            temp = struct('DrugName',MATLABsafename(ReplaceName({drugs_struct.name})), ...
                'SingleDoses',{drugs_struct.SingleDoses});
        end
        if isfield(drugs_struct,'Doses')
            for i=1:length(temp)
                temp(i).ComboDoses = drugs_struct(i).Doses;
            end
        end
        Drugs{iCL} = [Drugs{iCL} temp];
    end

    DrugDoses = zeros(height(t_CL{iCL}),length(DrugsNames{iCL}));
    for iDe=1:length(Designs)
        t_CL{iCL} = sortrows(t_CL{iCL},{'DesignSet' 'Replicate' 'Column' 'Row'});
        for iR = setdiff(unique(t_data.Replicate(strcmp(t_data.DesignSet, Designs{iDe}))), 0)'

            design_file = Barcode_table.DesignFile{strcmp(Barcode_table.CellLine,...
                char(CellLines(iCL))) & Barcode_table.Replicate==iR & strcmp(Barcode_table.DesignSet, Designs{iDe})};

            load([folder design_file ])
            temp_DrugNames = {drugs_struct.name};
            assert(all(ismember(temp_DrugNames, DrugsNames{iCL})))

            for iD=1:length(temp_DrugNames)
                DrugDoses(t_CL{iCL}.Replicate==iR & strcmp(t_CL{iCL}.DesignSet, Designs{iDe}),...
                    strcmp(temp_DrugNames{iD},DrugsNames{iCL})) = drugs_struct(iD).layout(:);
            end
        end
    end

    CtrlIdx = all(DrugDoses==0,2) & ~t_CL{iCL}.Day0;
    for iDe=1:length(Designs)
        for iR = setdiff(unique(t_CL{iCL}.Replicate), 0)'
            % remove corner as control to avoid bias (only the case if more
            % than 20 controls overall (enough control remain)
            CornerIdx = find(t_CL{iCL}.Replicate==iR & strcmp(t_CL{iCL}.DesignSet, Designs{iDe}) & ...
                ismember(t_CL{iCL}.Well,{'A1' 'A24' 'P1' 'P24'}));
            assert(length(CornerIdx)==4)
            if all(CtrlIdx(CornerIdx))
                CtrlIdx(CornerIdx)=false;
            end
        end
    end
    t_CL{iCL} = [t_CL{iCL} array2table([DrugDoses CtrlIdx], ...
        'VariableNames',{Drugs{iCL}.DrugName 'Ctrl'})];

    t_CL{iCL} = sortrows(t_CL{iCL},{'Replicate' 'Column' 'Untrt' ...
        'Day0' 'Ctrl'});
end

if savefile
    save([OutputFile '.mat'], 't_CL', 'Drugs', 'CellLines')
end
