function t_data = Read_ColumbusData(filename, barcode_file)
% t_data = Read_ColumbusData(filename, barcode_file)
%
%   import the results of plates scanned with Columbus
%
%   filename is the output of Columbus which should contain at least the headers:
%           Result (matching the barcodes)
%           Well
%           NumberOfAnalyzedFields
%           Nuclei_NumberOfObjects
%           URL (for checking the PlateID)
%
%   barcode_file should be a table with experimental design with headers:
%           Barcode
%           CellLine
%           Treatmentfile (refer to a table or .mat with with the
%               treatments on the plate - plate design) use   -   if untreated
%           DesignNumber:  replicate or treatment design
%           Time (in hours)
%           any Additional field with relevant plate properties to pass to the
%           data (collagen, ...)
%

t_raw = tsv2table(filename);


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

% correct well labels from A4 to A04
well_idx = cellfun(@length,t_raw.Well)==2;
t_raw.Well(well_idx) = strcat(cellcell2cellstr(regexp(t_raw.Well(well_idx),'^(\w)\d$','tokens')), ...
    '0',cellcell2cellstr(regexp(t_raw.Well(well_idx),'^\w(\d)$','tokens')));
%%

t_barcode = tsv2table(barcode_file);

assert(length(unique(t_raw.Barcode))>=height(t_barcode), ...
    'Found %i plates, but expected %i plates; missing %s', ...
    length(unique(t_raw.Barcode)), height(t_barcode), ...
    strjoin(setdiff(t_barcode.Barcode,unique(t_raw.Barcode))',' ') );

url = regexp(t_raw.URL,'/','split');
PlateID = vertcat(url{:});

PlateID = cellfun(@str2num,PlateID(:,end-1));

CellLine = cell(height(t_raw),1);
Barcode = cell(height(t_raw),1);
Treatmentfile = cell(height(t_raw),1);
DesignNumber = zeros(height(t_raw),1);
Time = zeros(height(t_raw),1);
Untrt = false(height(t_raw),1);
cnt = 0;
otherVariables = setdiff(varnames(t_barcode), {'Time' 'CellLine' 'Barcode' ...
            'Treatmentfile' 'DesignNumber' 'ExpNumber'});
for iBC = 1:height(t_barcode)

    idx = find(strcmp(t_raw.Barcode, t_barcode.Barcode{iBC}));
    assert(length(unique(PlateID(idx)))==1, '%i plate IDs found', ...
        length(unique(PlateID(idx))))
    assert(isempty(intersect(unique(PlateID(idx)), unique(PlateID(~idx)))))

    CellLine(idx) = t_barcode.CellLine(iBC);
    Barcode(idx) = t_barcode.Barcode(iBC);
    Treatmentfile(idx) = t_barcode.Treatmentfile(iBC);
    if ~isempty(strfind(t_barcode.Treatmentfile{iBC}, '.mat')) || ...
            isvariable(t_barcode,'DesignNumber')
        assert(isvariable(t_barcode,'DesignNumber'), ...
            'Treatment files .mat needs a DesignNumber column')
        DesignNumber(idx) = t_barcode.DesignNumber(iBC);
    else
        DesignNumber(idx) = NaN;
    end
    Untrt(idx) = strcmp(t_barcode.Treatmentfile(iBC),'-');
    Time(idx) = t_barcode.Time(iBC);

    for i = 1:length(otherVariables)
        eval([otherVariables{i} '(idx) = t_barcode.' otherVariables{i} '(iBC);'])
    end

    cnt = cnt+1;

end

for i = 1:length(otherVariables)
    eval([otherVariables{i} ' = ToColumn(' otherVariables{i} ');'])
    eval(['if length(' otherVariables{i} ')<height(t_raw), ' ...
        otherVariables{i} '(height(t_raw)+1) = ' otherVariables{i} '(1);' ...
        otherVariables{i} '(height(t_raw)+1) = []; end'])
end

assert(cnt<=length(unique(PlateID)))
if cnt<length(unique(PlateID))
    warning('some entries in the result file are unused!')
    Usedidx = ~cell2mat(cellfun2(@isempty,CellLine));
    Treatmentfile = Treatmentfile(Usedidx);
    CellLine = CellLine(Usedidx);
    DesignNumber = DesignNumber(Usedidx);
    Time = Time(Usedidx);
    for i = 1:length(otherVariables)
        eval([otherVariables{i} ' = ' otherVariables{i} '(Usedidx);'])
    end
else

    Usedidx = 1:height(t_raw);
end

Untrt = cellfun(@(x) strcmp(x,'-') || isempty(x), Treatmentfile);
assert(~any(DesignNumber==0 & ~Untrt), 'Some wells are not ''Untrt'' and don''t have a DesignNumber')
assert(all(Time(~Untrt)>0), 'Some treated wells don''t have a Time')

t_data = [table(CellLine,Treatmentfile,DesignNumber,Untrt,Time) ...
    t_raw(Usedidx, {'Well' 'Row' 'Column' 'Nuclei_NumberOfObjects' 'date'})];
if ~isempty(otherVariables)
    disp(['Added variables ''' cellstr2str(otherVariables, ''', ''') ''''])
    eval(['t_data = [t_data table(' cellstr2str(otherVariables, ',') ')];'])
end
t_data.Properties.VariableNames{'Nuclei_NumberOfObjects'} = 'Cellcount';
t_data = TableToCategorical(t_data,[1 2 6]);
