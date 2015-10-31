
function [mRNA_values, mRNA_labels] = ImportAnnotate_DGEdata(filename, plateinfo, Barcode, folder)
% [mRNA_values, mRNA_labels] = ImportAnnotate_DGEdata(filename, plateinfo, Barcode, folder)


p = struct('TimeCourse', false, 'Cellcount', 'none');
t_plateinfo = ImportCheckPlateInfo(plateinfo, p);


RNAfilename = [filename '.unq.refseq.umi.dat'];
rawRNA = importdata(RNAfilename);


Wells = CorrectWellsTo0xY(regexpcelltokens(rawRNA.textdata(1,2:end), '_([A-Q][0-9]*)'))';
t_labels = table(repmat({Barcode}, length(Wells),1), Wells, ...
    'variablenames', {'Barcode' 'Well'});

t_labels = AddPlateInfo_RawData(t_labels, t_plateinfo, '', p);
t_labels = Annotate_CellCountData(t_labels, folder);

INFOfilename = [filename '.unq.well_summary.dat'];
rawINFO = importdata(INFOfilename);
Wells = CorrectWellsTo0xY(regexpcelltokens(rawRNA.textdata(1,2:end), '_([A-Q][0-9]*)'))';
t_info = table(categorical(Wells), ...
    rawINFO.data(find(strcmp(rawINFO.textdata(:,1),'Refseq_Total'))-1,:)', ...
    rawINFO.data(find(strcmp(rawINFO.textdata(:,1),'Refseq_UMI'))-1,:)', ...
    'variablenames', {'Well' 'Refseq_Total' 'Refseq_UMI'});

t_labels = leftjoin(t_labels, t_info);

mRNA_names = table(categorical(rawRNA.textdata(2:end,1)),'VariableNames',{'gene'});
mRNA_values = rawRNA.data;

mRNA_labels = {mRNA_names t_labels};
